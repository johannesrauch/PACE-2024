#include "components.hpp"

#include <filesystem>

#include "bipartite_graph.hpp"
#include "input.hpp"

namespace fs = std::filesystem;

using T = uint16_t;

int main() {
    pace2024::bipartite_graph<T> graph;

    const fs::path filepath_matching = "tiny_test_set/matching_4_4.gr";
    pace2024::parse_input(filepath_matching, graph);
    {
        pace2024::components components(graph);
        assert(components.get_nof_components() == 4);
        components.build();
        for (std::size_t c = 0; c < components.get_nof_components(); ++c) {
            const pace2024::bipartite_graph<T>& component = components.get_component(c);
            assert(component.get_n_fixed() == 1);
            assert(component.get_n_free() == 1);
            assert(component.get_m() == 1);
        }
    }

    const fs::path filepath_star = "tiny_test_set/star_6.gr";
    pace2024::parse_input(filepath_star, graph);
    {
        pace2024::components components(graph);
        assert(components.get_nof_components() == 2);
        components.build();
        for (std::size_t c = 0; c < components.get_nof_components(); ++c) {
            const pace2024::bipartite_graph<T>& component = components.get_component(c);
            assert(component.get_n_fixed() == 1);
            assert(component.get_n_free() == 3);
            assert(component.get_m() == 3);
        }
    }

    const fs::path filepath_tree = "tiny_test_set/tree_6_10.gr";
    pace2024::parse_input(filepath_tree, graph);
    {
        pace2024::components components(graph);
        assert(components.get_nof_components() == 1);
        const pace2024::bipartite_graph<T>& component = components.get_component(0);
        assert(component.get_n_fixed() == 6);
        assert(component.get_n_free() == 10);
        assert(component.get_m() == 15);
    }

    graph.clear();
    {
        pace2024::components components(graph);
        assert(components.get_nof_components() == 0);
    }

    std::cout << "TEST::PACE2024::COMPONENTS:\t\tOK" << std::endl;
    return 0;
}