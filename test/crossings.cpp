#include "crossings.hpp"

#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>

#include "bipartite_graph.hpp"
#include "input.hpp"
#include "matrix.hpp"
#include "printf.hpp"
#include "random.hpp"
#include "test_utils.hpp"

namespace fs = std::filesystem;

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
    for (T u = 0; u < n_free; ++u) {
        for (T v = u + 1; v < n_free; ++v) {
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
                              const pace2024::folded_matrix<R>& folded_matrix) {
    std::vector<uint16_t> ordering(graph.get_n_free());
    pace2024::test::shuffle(ordering);

    std::clock_t start = std::clock();
    const R ref = pace2024::number_of_crossings(folded_matrix, ordering);
    std::clock_t end = std::clock();
    const double t_f = pace2024::test::time_in_ms(start, end);

    start = std::clock();
    uint32_t test = pace2024::number_of_crossings(graph, ordering);
    end = std::clock();
    assert(ref == test);
    const double t_g = pace2024::test::time_in_ms(start, end);

    pace2024::sparse_matrix<T, R> sparse_matrix(graph);
    start = std::clock();
    test = pace2024::number_of_crossings(sparse_matrix, ordering);
    end = std::clock();
    assert(ref == test);
    const double t_s = pace2024::test::time_in_ms(start, end);

    const char fastest = t_f <= t_g ? (t_f <= t_s ? 'f' : 's') : 'g';
    fmt::printf("%11.3f%11.3f%11.3f%11c\n", t_f, t_s, t_g, fastest);
}

/**
 * @brief tests crossings.hpp
 */
int main() {
    fmt::printf("%11s%11s%11s%11s%11s\n", "instance", "folded", "sparse", "graph", "fastest");
    pace2024::test::print_line(56);
    for (const auto& file : std::filesystem::directory_iterator("medium_test_set/instances")) {
        if (!file.is_regular_file()) continue;

        pace2024::bipartite_graph graph;
        pace2024::parse_input(file.path(), graph);
        pace2024::folded_matrix matrix(graph);

        fmt::printf("%11s", file.path().filename());
        test_crossing_numbers_of(graph, matrix);
        test_number_of_crossings(graph, matrix);
    }
    pace2024::test::print_line(56);
    std::cout << "TEST::PACE2024::CROSSING_NUMBER:\t\tOK" << std::endl;
    return 0;
}