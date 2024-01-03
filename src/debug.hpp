#ifndef PACE2024_DEBUG_HPP
#define PACE2024_DEBUG_HPP

#include "printf.hpp"

#ifndef NDEBUG
#define PACE2024_DEBUG_PRINTF(...) \
    do {                            \
        fmt::printf(__VA_ARGS__);          \
    } while (false)
#else
#define PACE2024_DEBUG_PRINTF(...) \
    do {                            \
    } while (false)
#endif

#endif