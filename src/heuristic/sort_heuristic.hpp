#ifndef PACE_HEURISTIC_SORT_HEURISTIC_HPP
#define PACE_HEURISTIC_SORT_HEURISTIC_HPP

#include "heuristic/shift_heuristic.hpp"
#include "utils/vector_utils.hpp"

namespace pace {

class highs_lp;

class sort_heuristic : public instance_view {
    shift_heuristic shift_h;

   public:
    sort_heuristic(instance &instance_)
        : instance_view(instance_), shift_h(instance_) {}

    sort_heuristic(const sort_heuristic &other) = delete;
    sort_heuristic(sort_heuristic &&other) = delete;
    sort_heuristic &operator=(const sort_heuristic &other) = delete;
    sort_heuristic &operator=(sort_heuristic &&other) = delete;

    crossing_number_t operator()(highs_lp &lp, std::vector<vertex_t> &ordering);
};

};  // namespace pace

#endif
