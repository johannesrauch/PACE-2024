#ifndef PACE_HEURISTIC_RELAX_HEURISTIC_HPP
#define PACE_HEURISTIC_RELAX_HEURISTIC_HPP

#include "exact/highs_lp.hpp"

namespace pace {

class relax_heuristic : public instance_view {
    highs_lp &lp;

   public:
    relax_heuristic(instance &instance_, highs_lp &lp) : instance_view(instance_), lp(lp) {}

    crossing_number_t operator()(std::vector<vertex_t> &ordering) {
        
    }
};

};  // namespace pace

#endif
