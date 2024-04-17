#ifndef PACE_LOG_DEBUG_PRINTF_HPP
#define PACE_LOG_DEBUG_PRINTF_HPP

#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>

#include "exact/info_structs.hpp"
#include "fmt/printf.hpp"

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
#define PACE_DEBUG_INFO(info_br, info_lp)                                                                           \
    do {                                                                                                            \
        if (info_br.n_iterations % 50 == 0) {                                                                       \
            fmt::printf("%11s%11s|%11s%11s%11s%11s|%11s%11s%11s%11s|%11s%11s%11s%11s\n", "lb", "ub", "iter",        \
                        "obj val", "nof rows", "nof nodes", "nof bucket", "nof added", "nof delete", "nof spared",  \
                        "uvw iter", "u", "v", "w");                                                                 \
        }                                                                                                           \
        fmt::printf("%11u%11u|%11u%11.0f%11u%11u|%11u%11u%11u%11u|%11u%11u%11u%11u\n", info_br.lower_bound,         \
                    info_br.upper_bound, info_br.n_iterations, info_lp.objective_value, info_lp.n_rows,             \
                    info_br.n_branch_nodes, info_lp.n_bucket_entries, info_lp.n_added_rows, info_lp.n_deleted_rows, \
                    info_lp.n_delete_rows_spared, info_lp.n_iterations_3cycles, info_lp.u_old, info_lp.v_old,       \
                    info_lp.w_old);                                                                                 \
        std::cout << std::flush;                                                                                    \
    } while (false)
#else
#define PACE_DEBUG_INFO(...) \
    do {                     \
    } while (false)
#endif

namespace pace {

namespace test {

class pipe_cout_to_file {
    std::ofstream out_f;
    std::streambuf* buf_old;

   public:
    pipe_cout_to_file(const char* filename) : out_f(filename), buf_old(std::cout.rdbuf(out_f.rdbuf())) {}
    ~pipe_cout_to_file() { std::cout.rdbuf(buf_old); }
};

};  // namespace test

};  // namespace pace

#endif