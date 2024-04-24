#ifndef PACE_LOG_DEBUG_PRINTF_HPP
#define PACE_LOG_DEBUG_PRINTF_HPP

#include "exact/info_structs.hpp"
#include "fmt/printf.hpp"

namespace pace {

namespace test {

void printf_info(const branch_and_cut_info& info_br, const highs_wrapper_info& info_lp) {
    if (info_br.n_iterations % 50 == 0) {
        fmt::printf("%11s%11s|%11s%11s%11s%11s%11s|%11s%11s%11s%11s%11s%11s|%11s%11s%11s%11s\n",  //
                    "lb", "ub",                                                                   //
                    "iter", "obj val", "n simplex", "n rows", "n nodes",                          //
                    "n bucket", "n added", "n delete", "n d slack", "n d luck", "n spared",       //
                    "uvw iter", "u", "v", "w");
    }
    fmt::printf("%11u%11u|%11u%11.1f%11u%11u%11u|%11u%11u%11u%11u%11u%11u|%11u%11u%11u%11u\n",  //
                info_br.lower_bound, info_br.upper_bound,                               //
                info_br.n_iterations, info_lp.objective_value, info_lp.n_iterations_simplex, info_lp.n_rows,
                info_br.n_branch_nodes,  //
                info_lp.n_bucket_entries, info_lp.n_added_rows, info_lp.n_deleted_rows, info_lp.n_deleted_rows_slack,
                info_lp.n_deleted_rows_bad_luck, info_lp.n_delete_rows_spared,  //
                info_lp.n_iterations_3cycles, info_lp.u_old, info_lp.v_old, info_lp.w_old);
}

};  // namespace test

};  // namespace pace

#ifdef PACE_DEBUG_PRINT
#define PACE_DEBUG_PRINTF(...)    \
    do {                          \
        fmt::printf(__VA_ARGS__); \
        std::cout << std::flush;  \
    } while (false)
#else
#define PACE_DEBUG_PRINTF(...) \
    do {                       \
    } while (false)
#endif

#ifdef PACE_DEBUG_PRINT
#define PACE_DEBUG_PRINTF_INFO(info_br, info_lp)          \
    do {                                           \
        pace::test::printf_info(info_br, info_lp); \
        std::cout << std::flush;                   \
    } while (false)
#else
#define PACE_DEBUG_PRINTF_INFO(...) \
    do {                     \
    } while (false)
#endif

#endif
