#include <cassert>
#include <filesystem>

#include "crossing_number.hpp"
#include "bipartite_graph.hpp"
#include "instance.hpp"
#include "matrix.hpp"

void testcase(std::string filepath) {
    pace2024::uint16_bipartite_graph graph(filepath);
    pace2024::folded_square_matrix<uint16_t> matrix(graph);

    for (uint16_t u = 0; u < graph.get_n1(); ++u) {
        for (uint16_t v = u + 1; v < graph.get_n1(); ++v) {
            auto [c_uv, c_vu] = pace2024::crossing_numbers_of<uint16_t, uint16_t>(graph, u, v);
            assert(matrix(u, v) == c_uv);
            assert(matrix(v, u) == c_vu);
        }
    }
}

int main() {
    for (const auto& file: std::filesystem::directory_iterator("tiny_test_set")) {
        testcase(file.path());
    }
    return 0;
}