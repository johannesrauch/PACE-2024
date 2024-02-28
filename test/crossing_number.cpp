#include "crossing_number.hpp"

#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>

#include "bipartite_graph.hpp"
#include "input.hpp"
#include "matrix.hpp"
#include "printf.hpp"
#include "random.hpp"

/**
 * @brief tests crossing_numbers_of function against crossing number matrix
 *
 * @tparam T vertex type
 * @tparam R return type
 * @param graph
 * @param matrix
 */
template <typename T, typename R>
void test_crossing_numbers_of(const pace2024::general_bipartite_graph<T>& graph,
                              const pace2024::folded_matrix<R>& matrix) {
    std::size_t n1 = graph.get_n1();
    for (uint16_t u = 0; u < n1; ++u) {
        for (uint16_t v = u + 1; v < n1; ++v) {
            auto [c_uv, c_vu] = pace2024::crossing_numbers_of<T, R>(graph, u, v);
            assert(matrix(u, v) == c_uv);
            assert(matrix(v, u) == c_vu);
        }
    }
}

/**
 * @brief tests the two crossing_number_of functions against each other
 *
 * @tparam T vertex type
 * @tparam R return type
 * @param graph
 * @param matrix
 */
template <typename T, typename R>
void test_crossing_number_of(const pace2024::general_bipartite_graph<T>& graph,
                             const pace2024::folded_matrix<R>& matrix) {
    std::vector<uint16_t> ordering(graph.get_n1());
    pace2024::test::shuffle(ordering);
    R ref = pace2024::crossing_number_of<T, R>(matrix, ordering);
    R tst = pace2024::crossing_number_of<T, R>(graph, ordering);
    // fmt::printf("%s,%s\n", ref, tst);
    assert(ref == tst);
}

/**
 * @brief tests crossing_number.hpp
 *
 * @return int
 */
int main() {
    for (const auto& file : std::filesystem::directory_iterator("tiny_test_set")) {
        pace2024::uint16_bipartite_graph graph(static_cast<const std::string>(file.path()));
        pace2024::uint16_folded_matrix matrix(graph);

        test_crossing_numbers_of(graph, matrix);
        test_crossing_number_of(graph, matrix);
    }

    std::cout << "TEST::PACE2024::CROSSING_NUMBER:\t\tOK" << std::endl;
    return 0;
}