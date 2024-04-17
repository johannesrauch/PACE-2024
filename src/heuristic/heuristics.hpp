#ifndef PACE_HEURISTICS_HEURISTICS_HPP
#define PACE_HEURISTICS_HEURISTICS_HPP

#include "heuristic/barycenter_heuristic.hpp"
#include "heuristic/median_heuristic.hpp"
#include "utils/crossings_utils.hpp"

namespace pace {

crossing_number_t heuristics(instance &instance, std::vector<vertex_t> &ordering) {
    PACE_DEBUG_PRINTF("start heuristic\n");

    const crossing_number_t &lb = instance.get_lower_bound();
    crossing_number_t ub = barycenter_heuristic{instance}(ordering);
    assert(lb <= ub);
    if (lb == ub) return ub;

    std::vector<vertex_t> ordering_(ordering.size());
    crossing_number_t ub_ = median_heuristic{instance}(ordering_);
    if (ub_ < ub) {
        ub = ub_;
        std::swap(ordering, ordering_);
    }
    assert(lb <= ub);
    if (lb == ub) return ub;

    ub = probmedian_heuristic{instance}(ordering, ub);
    PACE_DEBUG_PRINTF("end   heuristic\n");
    return ub;
}

};  // namespace pace

#endif
