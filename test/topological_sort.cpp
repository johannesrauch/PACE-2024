#include "utils/topological_sort.hpp"

#include <iostream>

#include "model/digraph.hpp"

/**
 * @brief tests utils/topological_sort.hpp (and model/digraph.hpp a little bit)
 */
int main() {
    {
        pace::general_digraph<uint8_t> digraph(6);
        digraph.add_arc(5, 2);
        digraph.add_arc(5, 0);
        digraph.add_arc(4, 0);
        digraph.add_arc(4, 1);
        digraph.add_arc(2, 3);
        digraph.add_arc(3, 1);
        std::vector<uint8_t> ordering;
        assert(pace::topological_sort(digraph, ordering));

        // for (uint8_t v: ordering) {
        //     std::cout << (int)v << std::endl;
        // }
        std::vector<uint8_t> ref_ordering = {5, 4, 2, 3, 1, 0};
        assert(ordering == ref_ordering);
    }

    {
        pace::general_digraph<uint8_t> digraph(6);
        digraph.add_arc(0, 1);
        digraph.add_arc(1, 2);
        digraph.add_arc(2, 3);
        digraph.add_arc(3, 4);
        digraph.add_arc(4, 5);
        digraph.add_arc(5, 0);
        std::vector<uint8_t> ordering;
        assert(!pace::topological_sort(digraph, ordering));
    }
    return 0;
}
