#ifndef PACE_HEURISTIC_RINS_HEURISTIC_HPP
#define PACE_HEURISTIC_RINS_HEURISTIC_HPP

#include "exact/highs_lp.hpp"
#include "model/instance.hpp"

namespace pace {

class rins_heuristic : public instance_view {
   public:
    rins_heuristic(instance &instance_) : instance_view(instance_) {}

    crossing_number_t operator()(const highs_lp &lp, std::vector<vertex_t> &ordering) {
        
    }
};

};  // namespace pace

#endif
