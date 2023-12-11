
#include "bipartite_graph.hpp"

#include <cassert>
#include <fstream>

#include "instance.hpp"

std::uint64_t ref_n0 = 4;
std::uint64_t ref_n1 = 4;
std::uint64_t ref_m = 8;
std::vector<std::pair<uint64_t, uint64_t>> ref_edges = {
    {1, 5}, {2, 5}, {2, 6}, {4, 6}, {1, 7}, {3, 7}, {3, 8}, {4, 8}};

int main() {
    std::ifstream input("tiny_test_set/cycle_8_sorted.gr", std::ios::in);
    // input >> std::cout.rdbuf();

    pace2024::instance instance;
    instance.parse(input);
    pace2024::bipartite_graph graph(instance);

    assert(ref_n0 == graph.get_n0());
    assert(ref_n1 == graph.get_n1());
    assert(ref_m == graph.get_m());

    uint64_t i = 0;
    auto adjacency_lists = graph.get_adjacency_lists();
    for (uint64_t y = 0; y < ref_n1; ++y) {
        for (uint64_t x : adjacency_lists[y]) {
            auto [ref_x, ref_y] = ref_edges[i];
            // std::cout << ref_x << " " << x + 1 << std::endl;
            assert(ref_x == x + 1);
            assert(ref_y == y + ref_n0 + 1);
            ++i;
        }
    }

    input.close();
    std::cout << "TEST::PACE2024::INSTANCE: OKAY" << std::endl;
    return 0;
}
