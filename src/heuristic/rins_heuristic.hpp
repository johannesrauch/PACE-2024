#ifndef PACE_HEURISTIC_RINS_HEURISTIC_HPP
#define PACE_HEURISTIC_RINS_HEURISTIC_HPP

#include "exact/highs_lp.hpp"
#include "heuristic/shift_heuristic.hpp"
#include "model/instance.hpp"

namespace pace {

class rins_heuristic : public instance_view {
    shift_heuristic shift_h;

   public:
    rins_heuristic(instance &instance_)
        : instance_view(instance_), shift_h(instance_) {}

    crossing_number_t operator()(highs_lp &lp, std::vector<vertex_t> &ordering);
};

};  // namespace pace

#endif
