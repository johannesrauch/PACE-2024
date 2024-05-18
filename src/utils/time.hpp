#ifndef PACE_UTILS_TLE_HPP
#define PACE_UTILS_TLE_HPP

#include <signal.h>

#include <chrono>

namespace pace {

volatile static sig_atomic_t tle = 0;

void time_limit_exceeded(int signum);

static std::chrono::time_point<std::chrono::system_clock> t0{
    std::chrono::system_clock::now()};

void reset_t0();

std::chrono::time_point<std::chrono::system_clock> now();

double elapsed_walltime_in_s(
    const std::chrono::time_point<std::chrono::system_clock> t1 = now(),
    const std::chrono::time_point<std::chrono::system_clock> t0 = t0);

};  // namespace pace

#endif
