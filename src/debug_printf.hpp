#ifndef PACE_DEBUG_PRINTF_HPP
#define PACE_DEBUG_PRINTF_HPP

#include <iostream>

#include "printf.hpp"

#ifdef NDEBUG_PRINT
#define PACE_DEBUG_PRINTF(...) \
    do {                           \
        fmt::printf(__VA_ARGS__);  \
        std::cout << std::flush;   \
    } while (false)
#else
#define PACE_DEBUG_PRINTF(...) \
    do {                           \
    } while (false)
#endif

#endif