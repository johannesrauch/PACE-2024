#ifndef PACE_LOG_DEBUG_PRINTF_HPP
#define PACE_LOG_DEBUG_PRINTF_HPP

#include "exact/info_structs.hpp"
#include "fmt/printf.hpp"

namespace pace {

namespace test {

void printf_info(const branch_and_cut_info& info_br, const highs_lp_info& info_lp);

void printf_lpinfo_line();

void printf_lpinfo(const highs_lp_info& info);

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

#endif
