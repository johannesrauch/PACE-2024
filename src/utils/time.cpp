#include "utils/time.hpp"

#include "log/debug_printf.hpp"

namespace pace {

timelimit &timelimit::singleton() {
    static timelimit timelimit;
    return timelimit;
}

std::chrono::time_point<std::chrono::system_clock> now() {
    return std::chrono::system_clock::now();
}

double elapsed_walltime_in_s(
    const std::chrono::time_point<std::chrono::system_clock> t1,
    const std::chrono::time_point<std::chrono::system_clock> t0) {
    return std::chrono::duration<double>(t1 - t0).count();
}

}  // namespace pace
