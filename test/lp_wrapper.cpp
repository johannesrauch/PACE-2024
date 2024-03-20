#include <filesystem>
#include <iostream>

#include "bipartite_graph.hpp"
#include "input.hpp"
#include "lp_wrapper.hpp"

int main() {
    const std::filesystem::path filepath = "tiny_test_set/instances/cycle_8_sorted.gr";
    pace::bipartite_graph graph;
    pace::parse_input(filepath, graph);
    pace::highs_wrapper lp(graph);
    std::cout << "TEST::PACE::LP_WRAPPER:\t\t\t\tOK" << std::endl;
}