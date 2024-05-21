#ifndef PACE_UTILS_TLE_HPP
#define PACE_UTILS_TLE_HPP

#include <chrono>
#include <csignal>

namespace pace {

/**
 * @brief singleton class
 */
class timelimit {
    timelimit() = default;
    static timelimit &singleton();
    std::sig_atomic_t tle{0};

    public:
    static bool was_sigterm_sent() {
        return singleton().tle;
    }

    static void sigterm_sent(int signum) {
        singleton().tle = 1;
    }
};

static std::chrono::time_point<std::chrono::system_clock> t0{
    std::chrono::system_clock::now()};

std::chrono::time_point<std::chrono::system_clock> now();

double elapsed_walltime_in_s(
    const std::chrono::time_point<std::chrono::system_clock> t1 = now(),
    const std::chrono::time_point<std::chrono::system_clock> t0 = t0);

};  // namespace pace

#endif
