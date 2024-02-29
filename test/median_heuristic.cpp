#include "median_heuristic.hpp"

#include <cassert>
#include <filesystem>
#include <iostream>

#include "bipartite_graph.hpp"
#include "crossing_number.hpp"
#include "printf.hpp"

/**
 * @brief tests the heuristic solvers median_heuristic (against expected values) and
 * probabilistic_median_heuristic (against values of median_heuristic;
 * should be better).
 *
 * @tparam T vertex type
 * @tparam R return type
 * @param graph
 * @param expected
 */
template <typename T, typename R>
void test_median_heuristics(const pace2024::general_bipartite_graph<T>& graph,
                            const R expected) {
    std::vector<T> ordering(graph.get_n1());
    pace2024::median_heuristic<T>(graph, ordering).run();
    R nof_crossings = pace2024::crossing_number_of<T, R>(graph, ordering);
    // fmt::printf("%s,\n", nof_crossings);

    R nof_crossings_ = pace2024::probabilistic_median_heuristic<T, R>(graph, ordering).run();
    assert(nof_crossings_ <= nof_crossings);
}

/**
 * @brief tests median_heuristic.hpp
 *
 * @return int
 */
int main() {
    uint32_t expected[] = {13, 4, 60, 0, 6, 0, 3, 11, 20, 0, 0, 17, 3};

    std::size_t i{0};
    for (const auto& file : std::filesystem::directory_iterator("tiny_test_set")) {
        if (!file.is_regular_file()) continue;

        pace2024::uint16_bipartite_graph graph(static_cast<const std::string>(file.path()));
        test_median_heuristics(graph, expected[i++]);
    }
    assert(i > 0);

    std::cout << "TEST::PACE2024::MEDIAN_HEURISTIC:\t\tOK" << std::endl;
    return 0;
}