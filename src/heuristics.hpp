#ifndef PACE_HEURISTICS_HPP
#define PACE_HEURISTICS_HPP

#include "crossings.hpp"
#include "median_heuristic.hpp"
#include "barycenter_heuristic.hpp"
#include "shift_heuristic.hpp"
#include "instance.hpp"

namespace pace {

template <typename T, typename R>
uint32_t heuristics(instance<T, R> &instance, std::vector<T> &ordering) {
    const uint32_t &lb = instance.get_lower_bound();
    uint32_t ub = barycenter_heuristic{instance}(ordering);
    assert(lb <= ub);
    if (lb == ub) return ub;

    std::vector<T> ordering_(ordering.size());
    uint32_t ub_ = median_heuristic{instance}(ordering_);
    if (ub_ < ub) {
        ub = ub_;
        std::swap(ordering, ordering_);
    }
    assert(lb <= ub);
    if (lb == ub) return ub;

    return probmedian_heuristic{instance}(ordering, ub);
}

};

#endif
