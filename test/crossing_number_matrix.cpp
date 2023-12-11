#include "crossing_number_matrix.hpp"

#include <iostream>

int main() {
    pace2024::instance instance("tiny_test_set/tree_6_10.gr");
    pace2024::bipartite_graph graph(instance);

    uint64_t n1 = graph.get_n1();
    pace2024::crossing_number_matrix ref_matrix(n1);
    pace2024::fill_naivly(graph, ref_matrix);
    pace2024::crossing_number_matrix matrix(graph);

    // pace2024::print_crossing_number_matrix(ref_matrix);
    // std::cout << std::endl;
    // pace2024::print_crossing_number_matrix(matrix);

    for (uint64_t i = 0; i < n1; ++i) {
        for (uint64_t j = 0; j < n1; ++j) {
            assert(ref_matrix(i, j) == matrix(i, j));
        }
    }

    std::cout << "TEST::PACE2024::CROSSING_NUMBER_MATRIX: OKAY" << std::endl;
}