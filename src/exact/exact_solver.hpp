#ifndef PACE_EXACT_EXACT_SOLVER_HPP
#define PACE_EXACT_EXACT_SOLVER_HPP

#include "exact/branch_and_cut.hpp"
#include "io/input.hpp"

namespace pace {

class exact_solver {
    input &in;

   public:
    exact_solver(input &in) : in(in) { in.try_split(); }

    crossing_number_t operator()(std::vector<vertex_t> &ordering);

   private:
    crossing_number_t solve_if_split(std::vector<vertex_t> &ordering);
};

}  // namespace pace

#endif
