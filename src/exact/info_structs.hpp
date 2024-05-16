#ifndef PACE_LOG_BRANCH_AND_CUT_INFO_HPP
#define PACE_LOG_BRANCH_AND_CUT_INFO_HPP

#include <chrono>

#include "types/types.hpp"

namespace pace {

static const std::chrono::time_point<std::chrono::system_clock> t0{
    std::chrono::system_clock::now()};

double elapsed_walltime_in_s(
    const std::chrono::time_point<std::chrono::system_clock> &t0 = t0);

struct highs_lp_info {
    const vertex_t &u_old;
    const vertex_t &v_old;
    const vertex_t &w_old;

    std::size_t n_cols{0};
    std::size_t n_rows{0};
    std::size_t n_runs{0};

    std::size_t n_bucket_entries{0};
    std::size_t n_added_rows{0};
    std::size_t n_deleted_rows{0};
    std::size_t n_delete_rows_spared{0};

    std::size_t n_solve_iters{0};
    std::size_t n_simplex_iters{0};
    std::size_t n_simplex_coldstart_iters{0};
    double n_avg_simplex_iters{0.};
    std::size_t n_3cycle_iters{0};
    std::size_t n_init_rows_candidates{0};
    bool new_3cycle_iter{false};

    double t_simplex{0.};
    double objective_value{0.};
    bool was_warmstart{false};

    double min_viol_score{0.};
    double max_viol_score{0.};
};

struct branch_and_cut_info {
    const crossing_number_t &lower_bound;
    const crossing_number_t &upper_bound;

    std::size_t n_iters{0};
    std::size_t n_search_nodes{0};
    std::size_t n_rows{0};
    std::size_t n_iter_3cycle_current_node{0};
    std::size_t depth;
    bool initially_solved{false};

    crossing_number_t n_crossings_h{0};

    double relax_h_confidence{-1.};
};

};  // namespace pace

#endif
