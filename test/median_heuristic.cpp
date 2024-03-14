#include "median_heuristic.hpp"

#include <cassert>
#include <filesystem>
#include <iostream>

#include "bipartite_graph.hpp"
#include "crossings.hpp"
#include "input.hpp"
#include "printf.hpp"

/**
 * @brief tests the heuristic solvers median_heuristic and
 * probabilistic_median_heuristic (against values of median_heuristic; should be better).
 */
template <typename T>
void test_median_heuristics(const pace2024::bipartite_graph<T>& graph) {
    std::vector<T> ordering(graph.get_n_free());
    pace2024::median_heuristic<T>(graph, ordering).run();
    uint32_t nof_crossings = pace2024::number_of_crossings(graph, ordering);
    uint32_t nof_crossings_ = pace2024::probabilistic_median_heuristic(graph, ordering).run();
    assert(nof_crossings_ <= nof_crossings);
}

/**
 * @brief tests median_heuristic.hpp
 */
int main() {
    for (const auto& file : std::filesystem::directory_iterator("tiny_test_set/instances")) {
        if (!file.is_regular_file()) continue;
        fmt::printf("%s\n", file.path());
        pace2024::bipartite_graph graph;
        pace2024::parse_input(file.path(), graph);
        test_median_heuristics(graph);
    }

    std::cout << "TEST::PACE2024::MEDIAN_HEURISTIC:\t\tOK" << std::endl;
    return 0;
}