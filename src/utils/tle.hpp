#include <signal.h>

#ifndef PACE_UTILS_TLE_HPP
#define PACE_UTILS_TLE_HPP

namespace pace {

static volatile sig_atomic_t tle = 0;

void time_limit_exceeded(int signum) {
    tle = 1;
}

};

#endif
