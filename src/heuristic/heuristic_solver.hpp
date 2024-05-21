#ifndef PACE_HEURISTIC_HEURISTIC_SOLVER_HPP
#define PACE_HEURISTIC_HEURISTIC_SOLVER_HPP

#include "io/input.hpp"
#include "heuristic/heuristics.hpp"

namespace pace {

class heuristic_solver {
    input &in;

    public:
    heuristic_solver(input &in) : in(in) {}

    crossing_number_t operator()(std::vector<vertex_t> &ordering) {
        instance instance_(in.get_graph());
        heuristics heuristic(instance_);
        return heuristic.uninformed(ordering, false);
    }
};

}  // namespace pace

#endif
