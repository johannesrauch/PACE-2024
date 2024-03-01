#ifndef PACE2024_DEBUG_PRINTF_HPP
#define PACE2024_DEBUG_PRINTF_HPP

#include <iostream>

#include "printf.hpp"

#ifdef NDEBUG_PRINT
#define PACE2024_DEBUG_PRINTF(...) \
    do {                           \
        fmt::printf(__VA_ARGS__);  \
        std::cout << std::flush;   \
    } while (false)
#else
#define PACE2024_DEBUG_PRINTF(...) \
    do {                           \
    } while (false)
#endif

#endif