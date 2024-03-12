#include <filesystem>
#include "bipartite_graph.hpp"
#include "input.hpp"
#include "components.hpp"
namespace fs = std::filesystem;

int main() {
    pace2024::bipartite_graph<uint16_t> graph;

    const fs::path filepath_matching = "tiny_test_set/matching_4_4.gr";
    pace2024::parse_input(filepath_matching, graph);
    {
        pace2024::components components(graph);
        assert(components.get_nof_components() == 4);
    }

    const fs::path filepath_star = "tiny_test_set/star_6.gr";
    pace2024::parse_input(filepath_star, graph);
    {
        pace2024::components components(graph);
        assert(components.get_nof_components() == 2);
    }

    const fs::path filepath_tree = "tiny_test_set/tree_6_10.gr";
    pace2024::parse_input(filepath_tree, graph);
    {
        pace2024::components components(graph);
        assert(components.get_nof_components() == 1);
    }

    graph.clear();
    {
        pace2024::components components(graph);
        assert(components.get_nof_components() == 0);
    }
    return 0;
}