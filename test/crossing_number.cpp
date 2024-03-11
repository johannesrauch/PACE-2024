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
void test_crossing_numbers_of(const pace2024::bipartite_graph<T>& graph,
                              const pace2024::folded_matrix<R>& matrix) {
    const std::size_t n_free = graph.get_n_free();
    for (uint16_t u = 0; u < n_free; ++u) {
        for (uint16_t v = u + 1; v < n_free; ++v) {
            auto [c_uv, c_vu] = pace2024::crossing_numbers_of<T, R>(graph, u, v);
            assert(matrix(u, v) == c_uv);
            assert(matrix(v, u) == c_vu);
        }
    }
}

/**
 * @brief tests the two number_of_crossings functions against each other
 *
 * @tparam T vertex type
 * @tparam R return type
 * @param graph
 * @param matrix
 */
template <typename T, typename R>
void test_number_of_crossings(const pace2024::bipartite_graph<T>& graph,
                             const pace2024::folded_matrix<R>& matrix) {
    std::vector<uint16_t> ordering(graph.get_n_free());
    pace2024::test::shuffle(ordering);
    R ref = pace2024::number_of_crossings<T, R>(matrix, ordering);
    R tst = pace2024::number_of_crossings<T, R>(graph, ordering);
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
        if (!file.is_regular_file()) continue;

        pace2024::uint16_bipartite_graph graph;
        pace2024::parse_input(file.path(), graph);
        pace2024::uint16_folded_matrix matrix(graph);

        test_crossing_numbers_of(graph, matrix);
        test_number_of_crossings(graph, matrix);
    }

    std::cout << "TEST::PACE2024::CROSSING_NUMBER:\t\tOK" << std::endl;
    return 0;
}