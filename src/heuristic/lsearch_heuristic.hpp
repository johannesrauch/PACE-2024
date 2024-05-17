#ifndef PACE_HEURISTIC_LSEARCH_HEURISTIC_HPP
#define PACE_HEURISTIC_LSEARCH_HEURISTIC_HPP

#include "model/instance.hpp"
#include "shift_heuristic.hpp"

namespace pace {

class lsearch_heuristic : public instance_view {
    uint16_t lsearch_width;
    shift_heuristic shift_h;

   public:
    lsearch_heuristic(instance &instance_, uint16_t lsearch_width = 16)
        : instance_view(instance_), lsearch_width(lsearch_width), shift_h(instance_) {}

    crossing_number_t operator()(std::vector<vertex_t> &ordering);
};

}  // namespace pace

#endif
