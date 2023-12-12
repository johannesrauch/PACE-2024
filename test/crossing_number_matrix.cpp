#include <iostream>

#include "matrix.hpp"

int main() {
    pace2024::instance instance("tiny_test_set/tree_6_10.gr");
    pace2024::bipartite_graph graph(instance);

    uint64_t n1 = graph.get_n1();
    pace2024::matrix ref_cr_matrix(n1, n1);
    pace2024::fill_naivly(graph, ref_cr_matrix);
    pace2024::matrix cr_matrix(graph);

    // pace2024::print_matrix(ref_cr_matrix);
    // std::cout << std::endl;
    // pace2024::print_matrix(cr_matrix);

    for (uint64_t i = 0; i < n1; ++i) {
        for (uint64_t j = 0; j < n1; ++j) {
            assert(ref_cr_matrix(i, j) == cr_matrix(i, j));
        }
    }

    std::cout << "TEST::PACE2024::CROSSING_NUMBER_MATRIX: OKAY" << std::endl;
    return 0;
}