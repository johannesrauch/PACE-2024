#ifndef PACE_HEURISTICS_HPP
#define PACE_HEURISTICS_HPP

#include "crossings.hpp"
#include "median_heuristic.hpp"
#include "barycenter_heuristic.hpp"
#include "shift_heuristic.hpp"
#include "instance.hpp"

namespace pace {

template <typename T, typename R>
uint32_t heuristics(const instance<T, R> &instance, std::vector<T> &ordering) {
    return heuristics(instance, ordering, get_lower_bound(instance.cr_matrix()));
}

template <typename T, typename R>
uint32_t heuristics(const instance<T, R> &instance, std::vector<T> &ordering, const uint32_t lb) {
    uint32_t ub = barycenter_heuristic{instance}(ordering);
    assert(lb <= ub);
    if (lb == ub) return ub;

    std::vector<T> ordering_(ordering.size());
    uint32_t ub_ = probmedian_heuristic{instance}(ordering_);
    if (ub_ < ub) {
        ub = ub_;
        ordering = ordering_;
    }

    assert(lb <= ub);
    return ub;
}

};

#endif
