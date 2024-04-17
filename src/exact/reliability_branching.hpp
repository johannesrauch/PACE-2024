#ifndef PACE_EXACT_RELIABILITY_BRANCHING_HPP
#define PACE_EXACT_RELIABILITY_BRANCHING_HPP

#ifndef PACE_CONST_RELIBRANCH_LOOKAHEAD
#define PACE_CONST_RELIBRANCH_LOOKAHEAD 8
#endif

#ifndef PACE_CONST_RELIBRANCH_RELIPARAM
#define PACE_CONST_RELIBRANCH_RELIPARAM 8
#endif

#include "lp_wrapper.hpp"

namespace pace {

class reliability_branching {
    highs_wrapper &lp;
    const std::size_t n_cols;
    std::vector<double> col_value;

    std::vector<double> up_cost_sums;
    std::vector<double> down_cost_sums;
    std::vector<crossing_number_t> n_up_summands;
    std::vector<crossing_number_t> n_down_summands;

    std::vector<double> scores;

   public:
    reliability_branching(highs_wrapper &lp)
        : lp(lp),  //
          n_cols(lp.get_nof_cols()),
          up_cost_sums(n_cols),
          down_cost_sums(n_cols),
          n_up_summands(n_cols),
          n_down_summands(n_cols),
          scores(n_cols) {}

    void update_up_cost(const std::size_t j, const double obj_val_new, const double obj_val_old, const double x_j) {
        assert(obj_val_new >= obj_val_old);
        assert(0. < x_j && x_j < 1.);
        up_cost_sums[j] += (obj_val_new - obj_val_old) / (1. - x_j);
        ++n_up_summands[j];
    }

    void update_down_cost(const std::size_t j, const double obj_val_new, const double obj_val_old, const double x_j) {
        assert(obj_val_new >= obj_val_old);
        assert(0. < x_j && x_j < 1.);
        down_cost_sums[j] += (obj_val_new - obj_val_old) / x_j;
        ++n_down_summands[j];
    }

    inline double score(const double &up_cost, const double &down_cost) {
        return 0.1667 * std::max(up_cost, down_cost) + 0.8333 * std::min(up_cost, down_cost);
    }

    void recompute_scores() {
        // scores of uninitialized columns is the average, which we compute first
        double up_cost_avg = 0.;
        double down_cost_avg = 0.;
        crossing_number_t up_cost_summands = 0;
        crossing_number_t down_cost_summands = 0;
        for (std::size_t j = 0; j < n_cols; ++j) {
            if (n_up_summands[j] > 0) {
                up_cost_avg += up_cost_sums[j];
                up_cost_summands += n_up_summands[j];
            }
            if (n_down_summands[j] > 0) {
                down_cost_avg += down_cost_sums[j];
                down_cost_summands += n_down_summands[j];
            }
        }
        if (up_cost_summands > 0) {
            up_cost_avg /= up_cost_summands;
        }
        if (down_cost_summands > 0) {
            down_cost_avg /= down_cost_summands;
        }

        // store scores
        for (std::size_t j = 0; j < n_cols; ++j) {
            const double x_j = col_value[j];
            double up_cost =             //
                (n_up_summands[j] > 0  //
                     ? up_cost_sums[j] / n_up_summands[j]
                     : up_cost_avg)  //
                * (1. - x_j);
            double down_cost =             //
                (n_down_summands[j] > 0  //
                     ? down_cost_sums[j] / n_down_summands[j]
                     : down_cost_avg)  //
                * x_j;
            scores[j] = score(up_cost, down_cost);
        }
    }

    double kinda_strong_branching(const std::size_t j) {
        const double obj_val_old = lp.get_objective_value();
        const double x_j = col_value[j];

        lp.fix_column(j, 1.);
        lp.run();
        assert(lp.is_optimal());
        double obj_val_new = lp.get_objective_value();
        const double delta_up = obj_val_new - obj_val_old;
        update_up_cost(j, obj_val_new, obj_val_old, x_j);
        
        lp.fix_column(j, 0.);
        lp.run();
        assert(lp.is_optimal());
        obj_val_new = lp.get_objective_value();
        const double delta_down = obj_val_new - obj_val_old;
        update_down_cost(j, obj_val_new, obj_val_old, x_j);

        lp.unfix_column(j);
        return score(delta_up, delta_down);
    }

    std::size_t get_branching_column() {
        // copy current column values
        lp.get_columns(col_value);

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

        // do strong branching for unreliable pseudo costs
        std::size_t best = *candidates.rbegin();
        double best_score = scores[best];
        std::size_t it_since_improve = 0;
        for (const std::size_t &j : candidates) {
            if (std::min(n_up_summands[j], n_down_summands[j]) < PACE_CONST_RELIBRANCH_RELIPARAM) {
                scores[j] = kinda_strong_branching(j);
                if (scores[j] > best_score) {
                    best = j;
                    best_score = scores[j];
                    it_since_improve = 0;
                } else {
                    ++it_since_improve;
                    if (it_since_improve >= PACE_CONST_RELIBRANCH_LOOKAHEAD) break;
                }
            }
        }

        return best;
    }

    double get_column_value(const std::size_t j) const {
        assert(j < n_cols);
        return col_value[j];
    }

    double get_up_cost(const std::size_t j) const {
        assert(j < n_cols);
        return up_cost_sums[j] / n_up_summands[j];
    }

    double get_down_cost(const std::size_t j) const {
        assert(j < n_cols);
        return down_cost_sums[j] / n_down_summands[j];
    }
};

};  // namespace pace

#endif
