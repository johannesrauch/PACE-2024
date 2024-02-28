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
        pace2024::general_digraph<uint8_t> digraph(6);
        digraph.add_arc(5, 2);
        digraph.add_arc(5, 0);
        digraph.add_arc(4, 0);
        digraph.add_arc(4, 1);
        digraph.add_arc(2, 3);
        digraph.add_arc(3, 1);
        std::vector<uint8_t> ordering;
        assert(pace2024::topological_sort(digraph, ordering));

        // for (uint8_t v: ordering) {
        //     std::cout << (int)v << std::endl;
        // }
        std::vector<uint8_t> ref_ordering = {5, 4, 2, 3, 1, 0};
        assert(ordering == ref_ordering);
    }

    {
        pace2024::general_digraph<uint8_t> digraph(6);
        digraph.add_arc(0, 1);
        digraph.add_arc(1, 2);
        digraph.add_arc(2, 3);
        digraph.add_arc(3, 4);
        digraph.add_arc(4, 5);
        digraph.add_arc(5, 0);
        std::vector<uint8_t> ordering;
        assert(!pace2024::topological_sort(digraph, ordering));
    }

    std::cout << "TEST::PACE2024::TOPOLOGICAL_SORT:\t\tOK" << std::endl;
    return 0;
}
