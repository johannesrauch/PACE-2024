#ifndef PACE_HEURISTICS_HPP
#define PACE_HEURISTICS_HPP

#include "median_heuristic.hpp"
#include "barycenter_heuristic.hpp"
#include "shift_heuristic.hpp"
#include "instance.hpp"

namespace pace {

template <typename T, typename R>
uint32_t heuristics(const instance<T, R> &instance, std::vector<T> &ordering) {
    uint32_t ub = probmedian_heuristic{instance}(ordering);

    std::vector<T> ordering_(ordering.size());
    uint32_t ub_ = barycenter_heuristic{instance}(ordering_);
    if (ub_ < ub) {
        ub = ub_;
        ordering = ordering_;
    }

    return ub;
}

};

#endif
