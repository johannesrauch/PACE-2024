#ifndef PACE_EXACT_RELIABILITY_BRANCHING_HPP
#define PACE_EXACT_RELIABILITY_BRANCHING_HPP

#include "exact/highs_lp.hpp"

namespace pace {

struct reliability_branching_params {
    const int32_t min_simplex_it{5000};
    const int32_t max_simplex_it{10000};
    const double max_simplex_it_factor{2};
    const uint8_t max_lookaheads{4};
    const uint8_t limit_reliable{4};
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
    reliability_branching(  //
        highs_lp &lp,
        reliability_branching_params params = reliability_branching_params())
        : lp(lp),  //
          n_cols(lp.get_n_cols()),
          up_cost_sums(n_cols),
          down_cost_sums(n_cols),
          n_up_samples(n_cols),
          n_down_samples(n_cols),
          scores(n_cols),
          params(params) {}

    void update_up_cost(const std::size_t j, const double obj_val_new,
                        const double obj_val_old, const double x_j);

    void update_down_cost(const std::size_t j, const double obj_val_new,
                          const double obj_val_old, const double x_j);

    inline double score(const double &up_cost, const double &down_cost) {
        return 0.1667 * std::max(up_cost, down_cost) +
               0.8333 * std::min(up_cost, down_cost);
    }

    void recompute_scores();

    double score_strong_branching(const std::size_t j);

    std::size_t get_branching_column();

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
