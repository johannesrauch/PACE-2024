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

    const bool first_graph_empty = in.is_first_graph_empty();
    if (first_graph_empty) {
        const std::size_t n_free = in.get_subgraph(0).get_n_free();
        identity(n_free, subordering);
        in.lift_ordering(0, subordering, ordering);
    }

    for (std::size_t i = first_graph_empty; i < n_subgraphs; ++i) {
        const bipartite_graph &subgraph{in.get_subgraph(i)};
        const std::size_t n_free = subgraph.get_n_free();
        assert(n_free > 0);
        if (n_free >= 2) {
            instance instance_(subgraph);
            branch_and_cut solver(instance_);
            n_crossings += solver(subordering);
        } else {
            identity(n_free, subordering);
        }
        in.lift_ordering(i, subordering, ordering);
    }

    assert(n_crossings == number_of_crossings(in.get_graph(), ordering));
    return n_crossings;
}

}  // namespace pace
