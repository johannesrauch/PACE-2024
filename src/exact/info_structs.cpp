#include "exact/info_structs.hpp"

namespace pace {

double elapsed_walltime_in_s(
    const std::chrono::time_point<std::chrono::system_clock> &t0) {
    const std::chrono::time_point<std::chrono::system_clock> now =
        std::chrono::system_clock::now();
    return std::chrono::duration<double>(now - t0).count();
}

};  // namespace pace
