#include "instance.hpp"

#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

std::uint64_t ref_n0 = 2;
std::uint64_t ref_n1 = 6;
std::uint64_t ref_m = 6;
std::vector<std::pair<uint64_t, uint64_t>> ref_edges = {{1, 3}, {2, 4}, {1, 5},
                                                        {2, 6}, {1, 7}, {2, 8}};

int main() {
    pace2024::instance instance("tiny_test_set/star_6.gr");

    assert(ref_n0 == instance.get_n0());
    assert(ref_n1 == instance.get_n1());
    assert(ref_m == instance.get_m());
    auto edges = instance.get_edges();
    for (uint64_t i = 0; i < ref_m; ++i) {
        assert(ref_edges[i].first == edges[i].first);
        assert(ref_edges[i].second == edges[i].second);
    }

    std::cout << "TEST::PACE2024::INSTANCE: OKAY" << std::endl;
    return 0;
}