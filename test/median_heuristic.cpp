#include <cassert>
#include <iostream>

#include "crossings.hpp"
#include "median_heuristic.hpp"
#include "instance.hpp"
#include "output.hpp"

int main() {
    pace2024::uint32_instance instance("tiny_test_set/website_20.gr");
    pace2024::uint32_bipartite_graph graph(instance);
    pace2024::folded_square_matrix<uint32_t> cr_matrix(graph);

    std::vector<uint32_t> ordering;
    pace2024::median_heuristic(graph, ordering);
    uint32_t nof_crossings = pace2024::compute_crossings(cr_matrix, ordering);
    assert(nof_crossings == 17);

    nof_crossings = pace2024::prob_median_heuristic(graph, cr_matrix, ordering);
    // std::cout << nof_crossings << std::endl;
    // pace2024::print_output(graph.get_n0(), ordering);

    std::cout << "TEST::PACE2024::MEDIAN_HEURISTIC: OKAY" << std::endl;
    return 0;
}