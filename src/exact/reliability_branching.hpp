#ifndef PACE_EXACT_RELIABILITY_BRANCHING_HPP
#define PACE_EXACT_RELIABILITY_BRANCHING_HPP

#include "exact/highs_lp.hpp"

namespace pace {

struct reliability_branching_params {
    const int32_t min_simplex_it{5000};
    const int32_t max_simplex_it{10000};
    const double max_simplex_it_factor{2};
    const uint8_t max_lookaheads{8};
    const uint8_t limit_reliable{8};
};

class reliability_branching {
    highs_lp &lp;
    const std::size_t n_cols;
    std::vector<double> col_value;
    double obj_val{0.};

    std::vector<double> up_cost_sums;
    std::vector<double> down_cost_sums;
    std::vector<uint32_t> n_up_samples;
    std::vector<uint32_t> n_down_samples;
    std::vector<double> scores;

    reliability_branching_params params;
    bool all_reliable{false};
    int32_t limit_simplex_it{0};

   public:
    reliability_branching(highs_lp &lp, reliability_branching_params params = reliability_branching_params())
        : lp(lp),  //
          n_cols(lp.get_n_cols()),
          up_cost_sums(n_cols),
          down_cost_sums(n_cols),
          n_up_samples(n_cols),
          n_down_samples(n_cols),
          scores(n_cols),
          params(params) {}

    void update_up_cost(const std::size_t j, const double obj_val_new, const double obj_val_old, const double x_j) {
        assert(j < n_cols);
        // assert(obj_val_new >= obj_val_old);
        up_cost_sums[j] += (obj_val_new - obj_val_old) / (1. - x_j);
        ++n_up_samples[j];
    }

    void update_down_cost(const std::size_t j, const double obj_val_new, const double obj_val_old, const double x_j) {
        assert(j < n_cols);
        // assert(obj_val_new >= obj_val_old);
        down_cost_sums[j] += (obj_val_new - obj_val_old) / x_j;
        ++n_down_samples[j];
    }

    inline double score(const double &up_cost, const double &down_cost) {
        return 0.1667 * std::max(up_cost, down_cost) + 0.8333 * std::min(up_cost, down_cost);
    }

    void recompute_scores() {
        // scores of uninitialized columns is the average, which we compute first
        double up_cost_avg = 0.;
        double down_cost_avg = 0.;
        crossing_number_t n_up_cost_summands = 0;
        crossing_number_t n_down_cost_summands = 0;

        for (std::size_t j = 0; j < n_cols; ++j) {
            if (n_up_samples[j] > 0) {
                up_cost_avg += up_cost_sums[j];
                n_up_cost_summands += n_up_samples[j];
            }
            if (n_down_samples[j] > 0) {
                down_cost_avg += down_cost_sums[j];
                n_down_cost_summands += n_down_samples[j];
            }
        }

        if (n_up_cost_summands > 0) {
            up_cost_avg /= n_up_cost_summands;
        }
        if (n_down_cost_summands > 0) {
            down_cost_avg /= n_down_cost_summands;
        }

        // store scores
        for (std::size_t j = 0; j < n_cols; ++j) {
            const double x_j = col_value[j];
            double up_cost =          //
                (n_up_samples[j] > 0  //
                     ? up_cost_sums[j] / n_up_samples[j]
                     : up_cost_avg)  //
                * (1. - x_j);
            double down_cost =          //
                (n_down_samples[j] > 0  //
                     ? down_cost_sums[j] / n_down_samples[j]
                     : down_cost_avg)  //
                * x_j;
            scores[j] = score(up_cost, down_cost);
        }
    }

    double score_strong_branching(const std::size_t j) {
        const double x_j = col_value[j];

        lp.fix_column(j, 1.);
        lp.run(limit_simplex_it);
        assert(lp.is_feasible());
        bool both_feasible = lp.is_feasible();
        double obj_val_new = lp.get_objective_value();
        const double delta_up = obj_val_new - obj_val;
        if (lp.is_feasible()) {
            update_up_cost(j, obj_val_new, obj_val, x_j);
        }

        lp.fix_column(j, 0.);
        lp.run(limit_simplex_it);
        both_feasible &= lp.is_feasible();
        obj_val_new = lp.get_objective_value();
        const double delta_down = obj_val_new - obj_val;
        if (lp.is_feasible()) {
            update_down_cost(j, obj_val_new, obj_val, x_j);
        }

        lp.unfix_column(j);
        if (both_feasible) {
            return score(delta_up, delta_down);
        } else {
            return scores[j];
        }
    }

    std::size_t get_branching_column() {
        lp.copy_column_values(col_value);
        obj_val = lp.get_objective_value();
        const int32_t max_simplex_it_suggested =
            lround(lp.get_info().n_iterations_simplex_avg * params.max_simplex_it_factor);
        limit_simplex_it = std::min(std::max(params.min_simplex_it, max_simplex_it_suggested), params.max_simplex_it);

        // collect nonintegral columns in candidates
        std::vector<std::size_t> candidates;
        for (std::size_t j = 0; j < n_cols; ++j) {
            if (!lp.is_column_integral(j)) {
                candidates.emplace_back(j);
            }
        }
        assert(candidates.size() > 0);

        // sort candidates s.t. score is non-increasing
        recompute_scores();
        std::sort(candidates.begin(), candidates.end(),
                  [&](const std::size_t &i, const std::size_t &j) -> bool { return scores[i] < scores[j]; });

        std::size_t best = *candidates.rbegin();
        double best_score = scores[best];
        if (all_reliable) return best;

        // do strong branching for unreliable pseudo costs
        all_reliable = true;
        std::size_t n_it_since_improve = 0;
        HighsInt frozen_basis_id = lp.freeze_basis();
        for (const std::size_t &j : candidates) {
            if (std::min(n_up_samples[j], n_down_samples[j]) >= params.limit_reliable) continue;

            all_reliable = false;
            scores[j] = score_strong_branching(j);
            if (scores[j] > best_score) {
                best = j;
                best_score = scores[j];
                n_it_since_improve = 1;
            } else {
                ++n_it_since_improve;
                if (n_it_since_improve > params.max_lookaheads) break;
            }
        }
        lp.unfreeze_basis(frozen_basis_id);

        PACE_DEBUG_PRINTF("chose %u with score %f\n", best, best_score);
        return best;
    }

    inline double get_column_value(const std::size_t j) const {
        assert(j < n_cols);
        return col_value[j];
    }

    inline double get_up_cost(const std::size_t j) const {
        assert(j < n_cols);
        return up_cost_sums[j] / n_up_samples[j];
    }

    inline double get_down_cost(const std::size_t j) const {
        assert(j < n_cols);
        return down_cost_sums[j] / n_down_samples[j];
    }
};

};  // namespace pace

#endif
