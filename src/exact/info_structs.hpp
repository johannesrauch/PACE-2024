#ifndef PACE_LOG_BRANCH_AND_CUT_INFO_HPP
#define PACE_LOG_BRANCH_AND_CUT_INFO_HPP

#include "types/types.hpp"

namespace pace {

struct highs_lp_info {
    const vertex_t &u_old;
    const vertex_t &v_old;
    const vertex_t &w_old;

    std::size_t n_cols{0};
    std::size_t n_rows{0};
    std::size_t n_runs{0};

    std::size_t n_deleted_rows{0};
    std::size_t n_deleted_rows_bad_luck{0};
    std::size_t n_deleted_rows_slack{0};
    std::size_t n_delete_rows_spared{0};
    
    std::size_t n_added_rows{0};
    
    std::size_t n_iterations_simplex{0};
    double n_iterations_simplex_avg{0.};
    std::size_t n_iterations_3cycles{0};
    std::size_t n_bucket_entries{0};
    std::size_t n_init_rows_candidates{0};

    double t_simplex{0.};
    double objective_value{0.};
};

struct branch_and_cut_info {
    const crossing_number_t &lower_bound;
    const crossing_number_t &upper_bound;

    std::size_t n_iterations{0};
    std::size_t n_branch_nodes{0};
    std::size_t n_rows{0};

    crossing_number_t n_crossings_h{0};
};

};

#endif
