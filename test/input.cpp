#include "input.hpp"
#include "test_utils.hpp"

#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

std::uint16_t ref_n_fixed = 2;
std::uint16_t ref_n_free = 6;
std::uint16_t ref_m = 6;
std::vector<std::pair<uint16_t, uint16_t>> ref_edges = {{1, 3}, {2, 4}, {1, 5},
                                                        {2, 6}, {1, 7}, {2, 8}};

/**
 * @brief tests input.hpp
 * 
 * @return int 
 */
int main() {
    const std::filesystem::path filepath{"tiny_test_set/instances/star_6.gr"};
    pace2024::bipartite_graph graph;
    pace2024::parse_input(filepath, graph);

    assert(ref_n_fixed == graph.get_n_fixed());
    assert(ref_n_free == graph.get_n_free());
    assert(ref_m == graph.get_m());
    const auto& edges = graph.get_edges();
    for (uint64_t i = 0; i < ref_m; ++i) {
        assert(ref_edges[i].first == edges[i].first + 1);
        assert(ref_edges[i].second == edges[i].second + 1 + ref_n_fixed);
    }

    std::cout << "TEST::PACE2024::INPUT:\t\t\t\tOK" << std::endl;
    return 0;
}