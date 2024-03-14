
#include "bipartite_graph.hpp"

#include <cassert>
#include <fstream>

#include "input.hpp"

namespace fs = std::filesystem;

constexpr uint64_t ref_n_fixed = 4;
constexpr uint64_t ref_n_free = 4;
constexpr uint64_t ref_m = 8;
const std::vector<std::pair<uint64_t, uint64_t>> ref_edges =
    {{1, 5}, {2, 5}, {2, 6}, {4, 6}, {1, 7}, {3, 7}, {3, 8}, {4, 8}};

/**
 * @brief tests bipartite_graph.hpp a little
 */
int main() {
    const fs::path filepath = "tiny_test_set/instances/cycle_8_sorted.gr";
    pace2024::uint64_bipartite_graph graph;
    pace2024::parse_input(filepath, graph);

    assert(ref_n_fixed == graph.get_n_fixed());
    assert(ref_n_free == graph.get_n_free());
    assert(ref_m == graph.get_m());

    uint64_t i = 0;
    for (uint64_t y = 0; y < ref_n_free; ++y) {
        for (const uint64_t& x : graph.get_neighbors_of_free(y)) {
            auto [ref_x, ref_y] = ref_edges[i];
            // std::cout << ref_x << " " << x + 1 << std::endl;
            assert(ref_x == x + 1);
            assert(ref_y == y + ref_n_fixed + 1);
            ++i;
        }
    }

    std::cout << "TEST::PACE2024::BIPARTITE_GRAPH:\t\tOK" << std::endl;
    return 0;
}
