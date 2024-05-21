#include "exact/exact_solver.hpp"

namespace pace {

//
// public methods
//

crossing_number_t exact_solver::operator()(std::vector<vertex_t> &ordering) {
    PACE_DEBUG_PRINTF(  //
        "start exact_solver, n_subgraphs=%u\n", in.get_n_subgraphs());
    crossing_number_t n_cr;
    if (in.get_n_subgraphs() >= 2) {
        n_cr = solve_if_split(ordering);
    } else {
        instance instance_(in.get_graph());
        branch_and_cut solver(instance_);
        n_cr = solver(ordering);
    }
    PACE_DEBUG_PRINTF("end   exact_solver\n\n");
    return n_cr;
}

//
// private methods
//

crossing_number_t exact_solver::solve_if_split(
    std::vector<vertex_t> &ordering) {
    const std::size_t n_subgraphs = in.get_n_subgraphs();
    assert(n_subgraphs >= 2);
    std::vector<vertex_t> subordering;
    crossing_number_t n_crossings{0};
    ordering.resize(in.get_graph().get_n_free());

    for (std::size_t i = 0; i < n_subgraphs; ++i) {
        const bipartite_graph &subgraph{in.get_subgraph(i)};
        n_crossings += solve(subgraph, subordering);
        in.lift_ordering(i, subordering, ordering);
    }

    assert(n_crossings == number_of_crossings(in.get_graph(), ordering));
    return n_crossings;
}

crossing_number_t exact_solver::solve(  //
    const bipartite_graph &graph, std::vector<vertex_t> &ordering) {
    const std::size_t n_free = graph.get_n_free();
    
    if (graph.get_m() == 0 || n_free == 1) {
        identity(n_free, ordering);
        return 0;
    }

    instance instance_(graph);
    if (n_free == 2) {
        ordering.resize(2);
        const auto [c_01, c_10] = instance_.get_cr_numbers(0, 1);
        if (c_01 < c_10) {
            ordering[0] = 0;
            ordering[1] = 1;
            return c_01;
        } else {
            ordering[0] = 1;
            ordering[1] = 0;    
            return c_10;
        }
    }

    branch_and_cut solver(instance_);
    return solver(ordering);
}

}  // namespace pace
