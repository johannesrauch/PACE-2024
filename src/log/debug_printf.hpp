#ifndef PACE_LOG_DEBUG_PRINTF_HPP
#define PACE_LOG_DEBUG_PRINTF_HPP

#include <vector>

#include "exact/info_structs.hpp"
#include "fmt/printf.hpp"
#include "utils/time.hpp"

namespace pace {

namespace test {

void printf_info(const branch_and_cut_info& info_br,
                 const highs_lp_info& info_lp);

void printf_lpinfo_line();

void printf_lpinfo(const highs_lp_info& info);

void printf_bounds(const crossing_number_t& lb, const crossing_number_t& ub);

void printf_graph(  //
    const std::size_t n0, const std::size_t n1, const std::size_t m);

void printf_summary(const branch_and_cut_info& info);

template <typename T>
static void printf_vector(const std::vector<T>& vec) {
    for (const T& v : vec) {
        fmt::printf("%s ", v);
    }
    fmt::printf("\n");
}

template <typename T>
static void printf_vector(const std::vector<std::pair<T, T>>& vec) {
    for (const std::pair<T, T>& v : vec) {
        fmt::printf("(%s,%s) ", v.first, v.second);
    }
    fmt::printf("\n");
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
#define PACE_DEBUG_PRINTF_INFO(info_br, info_lp)   \
    do {                                           \
        pace::test::printf_info(info_br, info_lp); \
        std::cout << std::flush;                   \
    } while (false)
#else
#define PACE_DEBUG_PRINTF_INFO(...) \
    do {                            \
    } while (false)
#endif

#ifdef PACE_DEBUG_PRINT
#define PACE_DEBUG_PRINTF_LPINFO_LINE()   \
    do {                                  \
        pace::test::printf_lpinfo_line(); \
    } while (false)
#else
#define PACE_DEBUG_PRINTF_LPINFO_LINE() \
    do {                                \
    } while (false)
#endif

#ifdef PACE_DEBUG_PRINT
#define PACE_DEBUG_PRINTF_LPINFO(info)   \
    do {                                 \
        pace::test::printf_lpinfo(info); \
    } while (false)
#else
#define PACE_DEBUG_PRINTF_LPINFO(...) \
    do {                              \
    } while (false)
#endif

#ifdef PACE_DEBUG_PRINT
#define PACE_DEBUG_PRINTF_BOUNDS(lb, ub)   \
    do {                                   \
        pace::test::printf_bounds(lb, ub); \
    } while (false)
#else
#define PACE_DEBUG_PRINTF_BOUNDS(...) \
    do {                              \
    } while (false)
#endif

#ifdef PACE_DEBUG_PRINT
#define PACE_DEBUG_PRINTF_GRAPH(n0, n1, m)   \
    do {                                     \
        pace::test::printf_graph(n0, n1, m); \
    } while (false)
#else
#define PACE_DEBUG_PRINTF_GRAPH(...) \
    do {                             \
    } while (false)
#endif

#ifdef PACE_DEBUG_PRINT
#define PACE_DEBUG_PRINTF_SUMMARY(info)   \
    do {                                  \
        pace::test::printf_summary(info); \
    } while (false)
#else
#define PACE_DEBUG_PRINTF_SUMMARY(...) \
    do {                               \
    } while (false)
#endif

#ifdef PACE_DEBUG_PRINT
#define PACE_DEBUG_PRINTF_VECTOR(vec)   \
    do {                                \
        pace::test::printf_vector(vec); \
    } while (false)
#else
#define PACE_DEBUG_PRINTF_VECTOR(...) \
    do {                              \
    } while (false)
#endif

#endif
