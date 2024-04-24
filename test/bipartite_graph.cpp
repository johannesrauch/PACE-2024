
#include "io/parse_input.hpp"

namespace fs = std::filesystem;

constexpr uint64_t ref_n_fixed = 4;
constexpr uint64_t ref_n_free = 4;
constexpr uint64_t ref_m = 8;
const std::vector<std::pair<uint64_t, uint64_t>> ref_edges =
    {{1, 5}, {2, 5}, {2, 6}, {4, 6}, {1, 7}, {3, 7}, {3, 8}, {4, 8}};

/**
 * @brief tests model/bipartite_graph.hpp a little
 */
int main() {
    const fs::path filepath = "tiny_test_set/instances/cycle_8_sorted.gr";
    pace::bipartite_graph graph;
    pace::parse_input(filepath, graph);

    assert(ref_n_fixed == graph.get_n_fixed());
    assert(ref_n_free == graph.get_n_free());
    assert(ref_m == graph.get_m());

    using vertex_t = pace::bipartite_graph::vertex_t;
    vertex_t i = 0;
    for (vertex_t y = 0; y < ref_n_free; ++y) {
        for (const vertex_t& x : graph.get_neighbors(y)) {
            auto [ref_x, ref_y] = ref_edges[i];
            // std::cout << ref_x << " " << x + 1 << std::endl;
            assert(ref_x == x + 1);
            assert(ref_y == y + ref_n_fixed + 1);
            ++i;
        }
    }
    return 0;
}
