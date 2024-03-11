#include "input.hpp"

#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

std::uint64_t ref_n_fixed = 2;
std::uint64_t ref_n_free = 6;
std::uint64_t ref_m = 6;
std::vector<std::pair<uint64_t, uint64_t>> ref_edges = {{1, 3}, {2, 4}, {1, 5},
                                                        {2, 6}, {1, 7}, {2, 8}};

/**
 * @brief tests input.hpp
 * 
 * @return int 
 */
int main() {
    pace2024::uint64_input input("tiny_test_set/star_6.gr");

    assert(ref_n_fixed == input.get_n0());
    assert(ref_n_free == input.get_n1());
    assert(ref_m == input.get_m());
    auto edges = input.get_edges();
    for (uint64_t i = 0; i < ref_m; ++i) {
        assert(ref_edges[i].first == edges[i].first);
        assert(ref_edges[i].second == edges[i].second);
    }

    std::cout << "TEST::PACE2024::INPUT:\t\t\t\tOK" << std::endl;
    return 0;
}