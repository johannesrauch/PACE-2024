#include <iostream>

#include "matrix.hpp"

int main() {
    pace2024::uint64_instance instance("tiny_test_set/tree_6_10.gr");
    pace2024::uint64_bipartite_graph graph(instance);

    uint64_t n1 = graph.get_n1();
    pace2024::folded_square_matrix<uint64_t> ref_cr_matrix(n1);
    pace2024::compute_crossing_numbers_naivly(graph, ref_cr_matrix);
    pace2024::matrix<uint64_t> cr_matrix(graph);

    // pace2024::print_matrix(ref_cr_matrix);
    // std::cout << std::endl;
    // pace2024::print_matrix(cr_matrix);

    for (std::size_t i = 0; i < n1; ++i) {
        for (std::size_t j = 0; j < n1; ++j) {
            if (i == j) continue;
            assert(ref_cr_matrix(i, j) == cr_matrix(i, j));
        }
    }

    pace2024::folded_square_matrix<uint64_t> fl_matrix(graph);
    for (uint64_t i = 0; i < n1; ++i) {
        for (uint64_t j = 0; j < n1; ++j) {
            if (i == j) continue;
            // std::cout << ref_cr_matrix(i, j) << "," << fl_matrix(i, j)
            //         << std::endl;
            assert(ref_cr_matrix(i, j) == fl_matrix(i, j));
        }
    }

    std::cout << "TEST::PACE2024::CROSSING_NUMBER_MATRIX: OKAY" << std::endl;
    return 0;
}