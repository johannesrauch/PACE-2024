#include "lp_wrapper.hpp"

int main() {
    const std::string filepath = "tiny_test_set/cycle_8_sorted.gr";
    pace2024::uint16_bipartite_graph graph(filepath);
    pace2024::clp_wrapper lp(graph);
    std::cout << "TEST::PACE2024::LP_WRAPPER:\t\t\t\tOK" << std::endl;
}