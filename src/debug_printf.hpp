#ifndef PACE_DEBUG_PRINTF_HPP
#define PACE_DEBUG_PRINTF_HPP

#include <iostream>

#include "printf.hpp"

#ifdef DEBUG_PRINT
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

#ifdef DEBUG_PRINT
#define PACE_DEBUG_INFO(info_br, info_lp)                                                                       \
    do {                                                                                                        \
        if (info_br.nof_iterations % 50 == 0) {                                                                 \
            fmt::printf("%11s%11s|%11s%11s%11s|%11s%11s|%11s%11s%11s\n", "lb", "ub", "nof iter", "nof rows",    \
                        "nof nodes", "nof added r", "nof del r", "u", "v", "w");                                \
        }                                                                                                       \
        fmt::printf("%11u%11u|%11u%11u%11u|%11u%11u|%11u%11u%11u\n", info_br.lower_bound, info_br.upper_bound,  \
                    info_br.nof_iterations, info_lp.nof_rows, info_br.nof_branch_nodes, info_lp.nof_added_rows, \
                    info_lp.nof_deleted_rows, info_lp.u_old, info_lp.v_old, info_lp.w_old);                     \
        std::cout << std::flush;                                                                                \
    } while (false)
#else
#define PACE_DEBUG_INFO(...) \
    do {                     \
    } while (false)
#endif

#endif