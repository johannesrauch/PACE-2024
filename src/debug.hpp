#ifndef PACE2024_DEBUG_HPP
#define PACE2024_DEBUG_HPP

#include "printf.hpp"

#ifndef NDEBUG
#define PACE2024_DEBUG_PRINTF(args) \
    do {                            \
        fmt::printf(args);          \
    } while (false)
#else
#define PACE2024_DEBUG_PRINTF(args) \
    do {                            \
    } while (false)
#endif

#endif