#include "matrix.hpp"

#include <iostream>
#include <filesystem>

/**
 * @brief tests the different matrices and the different
 * ways to compute the crossing number matrix.
 * 
 * @tparam T vertex type
 * @param graph 
 */
template<typename T>
void test_matrix_and_compute_crossing_numbers(const pace2024::general_bipartite_graph<T>& graph) {
    const std::size_t n1{graph.get_n1()};
    pace2024::uint64_folded_matrix ref_matrix(n1);
    pace2024::uint64_matrix tst_matrix_1(graph);
    pace2024::uint64_folded_matrix tst_matrix_2(graph);
    pace2024::test::compute_crossing_numbers_naivly(graph, ref_matrix);
    
    assert(pace2024::test::equals(ref_matrix, tst_matrix_1, true));
    assert(pace2024::test::equals(ref_matrix, tst_matrix_2, true));
}

/**
 * @brief tests matrix.hpp; its classes and functions
 * 
 * @return int 
 */
int main() {
    for (const auto& file : std::filesystem::directory_iterator("tiny_test_set")) {
        pace2024::uint16_bipartite_graph graph(static_cast<const std::string>(file.path()));
        test_matrix_and_compute_crossing_numbers(graph);        
    }

    std::cout << "TEST::PACE2024::CROSSING_NUMBER_MATRIX: OKAY" << std::endl;
    return 0;
}