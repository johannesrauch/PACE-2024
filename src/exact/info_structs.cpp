#include "exact/info_structs.hpp"

namespace pace {

void reset_t0() { t0 = now(); }

std::chrono::time_point<std::chrono::system_clock> now() {
    return std::chrono::system_clock::now();
}

double elapsed_walltime_in_s(
    const std::chrono::time_point<std::chrono::system_clock> t1,
    const std::chrono::time_point<std::chrono::system_clock> t0) {
    return std::chrono::duration<double>(t1 - t0).count();
}

};  // namespace pace
