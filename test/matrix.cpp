#include "matrix.hpp"

#include <iostream>

#include "input.hpp"
#include "test_utils.hpp"

namespace fs = std::filesystem;

/**
 * @brief tests the different matrices and the different
 * ways to compute the crossing number matrix
 */
void test_crossing_matrices(const fs::path filepath) {
    pace::bipartite_graph graph;
    pace::parse_input(filepath, graph);
    const std::size_t n_free{graph.get_n_free()};

    // reference; naive computation
    pace::folded_matrix ref_matrix(n_free);
    std::clock_t start = std::clock();
    pace::test::fill_crossing_matrix_naivly(graph, ref_matrix);
    std::clock_t end = std::clock();
    const double t_r = pace::test::time_in_ms(start, end);

    // normal matrix
    start = std::clock();
    pace::matrix matrix(graph);
    end = std::clock();
    const double t_m = pace::test::time_in_ms(start, end);

    // cache optimized matrix
    start = std::clock();
    pace::folded_matrix folded_matrix(graph);
    end = std::clock();
    const double t_f = pace::test::time_in_ms(start, end);

    // sparse matrix
    start = std::clock();
    pace::sparse_matrix sparse_matrix(graph);
    end = std::clock();
    assert(sparse_matrix.good());
    const double t_s = pace::test::time_in_ms(start, end);

    // test equality
    assert(pace::test::equals(ref_matrix, matrix, true));
    assert(pace::test::equals(ref_matrix, folded_matrix, true));
    assert(pace::test::equals(ref_matrix, sparse_matrix, true));

    const char fastest = t_f <= t_m ? (t_f <= t_s ? 'f' : 's') : (t_m <= t_s ? 'm' : 's');
    fmt::printf("%11s%11.3f%11.3f%11.3f%11.3f%11c\n",  //
                filepath.filename(), t_r, t_m, t_f, t_s, fastest);
}

void test_crossing_matrices_with(const fs::path dirpath) {
    fmt::printf("start test_crossing_matrices_with(%s)\n\n", dirpath);
    fmt::printf("%11s%11s%11s%11s%11s%11s\n",  //
                "instance", "reference", "matrix", "folded", "sparse", "fastest");
    pace::test::print_line(67);
    for (const auto& file : std::filesystem::directory_iterator("medium_test_set/instances")) {
        if (!file.is_regular_file()) continue;
        test_crossing_matrices(file.path());
        std::cout << std::flush;
    }
    fmt::printf("\nend   test_crossing_matrices_with(%s)\n\n", dirpath);
}

/**
 * @brief tests matrix.hpp; its classes and functions
 *
 * @return int
 */
int main(int argc, char** argv) {
    fs::path dirpath;
    if (argc == 1) {
        dirpath = "medium_test_set/instances";
    } else {
        dirpath = "big_test_set/instances";
    }
    test_crossing_matrices_with(dirpath);

    std::cout << "TEST::PACE::MATRIX:\t\t\t\t\tOK" << std::endl;
    return 0;
}