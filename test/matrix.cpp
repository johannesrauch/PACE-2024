#include "matrix.hpp"

#include <iostream>

#include "input.hpp"
#include "test_utils.hpp"

namespace fs = std::filesystem;

/**
 * @brief tests the different matrices and the different
 * ways to compute the crossing number matrix
 */
void test_matrix_and_compute_crossing_numbers(const fs::path filepath) {
    pace2024::bipartite_graph graph;
    pace2024::parse_input(filepath, graph);
    const std::size_t n_free{graph.get_n_free()};

    // reference; naive computation
    pace2024::folded_matrix ref_matrix(n_free);
    std::clock_t start = std::clock();
    pace2024::test::fill_crossing_matrix_naivly(graph, ref_matrix);
    std::clock_t end = std::clock();
    const double t_r = pace2024::test::time_in_ms(start, end);

    // normal matrix
    start = std::clock();
    pace2024::matrix matrix(graph);
    end = std::clock();
    const double t_m = pace2024::test::time_in_ms(start, end);

    // cache optimized matrix
    start = std::clock();
    pace2024::folded_matrix folded_matrix(graph);
    end = std::clock();
    const double t_f = pace2024::test::time_in_ms(start, end);

    // sparse matrix
    start = std::clock();    
    pace2024::sparse_matrix sparse_matrix(graph);
    end = std::clock();
    assert(sparse_matrix.good());
    const double t_s = pace2024::test::time_in_ms(start, end);

    // test equality
    assert(pace2024::test::equals(ref_matrix, matrix, true));
    assert(pace2024::test::equals(ref_matrix, folded_matrix, true));
    assert(pace2024::test::equals(ref_matrix, sparse_matrix, true));

    const char fastest = t_f <= t_m ? (t_f <= t_s ? 'f' : 's') : 'm';
    fmt::printf("%11s%11.3f%11.3f%11.3f%11.3f%11c\n", filepath.filename(), t_r, t_m, t_f, t_s, fastest);
}

/**
 * @brief tests matrix.hpp; its classes and functions
 *
 * @return int
 */
int main() {
    fmt::printf("%11s%11s%11s%11s%11s%11s\n", "instance", "reference", "matrix", "folded", "sparse", "fastest");
    pace2024::test::print_line(67);
    for (const auto& file : std::filesystem::directory_iterator("medium_test_set/instances")) {
        if (!file.is_regular_file()) continue;
        test_matrix_and_compute_crossing_numbers(file.path());
    }
    pace2024::test::print_line(67);
    std::cout << "TEST::PACE2024::MATRIX:\t\t\t\t\tOK" << std::endl;
    return 0;
}