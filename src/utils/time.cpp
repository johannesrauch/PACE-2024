#include "utils/time.hpp"

namespace pace {

void time_limit_exceeded(int signum) { tle = 1; }

void reset_t0() { t0 = now(); }

std::chrono::time_point<std::chrono::system_clock> now() {
    return std::chrono::system_clock::now();
}

double elapsed_walltime_in_s(
    const std::chrono::time_point<std::chrono::system_clock> t1,
    const std::chrono::time_point<std::chrono::system_clock> t0) {
    return std::chrono::duration<double>(t1 - t0).count();
}

}  // namespace pace
