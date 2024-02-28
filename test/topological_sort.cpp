#include "topological_sort.hpp"

#include <cassert>
#include <iostream>

#include "graph.hpp"

/**
 * @brief tests topological_sort.hpp (and graph.hpp a little bit)
 *
 * @return int
 */
int main() {
    {
        pace2024::general_digraph<uint8_t> graph(6);
        graph.add_edge(5, 2);
        graph.add_edge(5, 0);
        graph.add_edge(4, 0);
        graph.add_edge(4, 1);
        graph.add_edge(2, 3);
        graph.add_edge(3, 1);
        std::vector<uint8_t> ordering;
        assert(pace2024::topological_sort(graph, ordering));

        // for (uint8_t v: ordering) {
        //     std::cout << (int)v << std::endl;
        // }
        std::vector<uint8_t> ref_ordering = {5, 4, 2, 3, 1, 0};
        assert(ordering == ref_ordering);
    }

    {
        pace2024::general_digraph<uint8_t> graph(6);
        graph.add_edge(0, 1);
        graph.add_edge(1, 2);
        graph.add_edge(2, 3);
        graph.add_edge(3, 4);
        graph.add_edge(4, 5);
        graph.add_edge(5, 0);
        std::vector<uint8_t> ordering;
        assert(!pace2024::topological_sort(graph, ordering));
    }

    std::cout << "TEST::PACE2024::TOPOLOGICAL_SORT: OKAY" << std::endl;
    return 0;
}
