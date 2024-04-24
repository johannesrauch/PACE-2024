#include "matrix/matrix.hpp"

#include "io/parse_input.hpp"
#include "utils/test_utils.hpp"

namespace fs = std::filesystem;

/**
 * @brief tests the different matrices and the different ways to compute the crossing number matrix
 */
void test_crossing_matrix(const fs::path filepath) {
    pace::bipartite_graph graph;
    pace::parse_input(filepath, graph);
    const std::size_t n_free{graph.get_n_free()};

    // reference; naive computation
    pace::crossing_matrix ref_matrix(n_free);
    std::clock_t start = std::clock();
    pace::test::fill_crossing_matrix_naivly(graph, ref_matrix);
    std::clock_t end = std::clock();
    const double t_r = pace::test::time_in_ms(start, end);

    // normal matrix
    start = std::clock();
    pace::matrix<uint32_t> matrix(graph);
    end = std::clock();
    const double t_m = pace::test::time_in_ms(start, end);

    // cache optimized matrix
    start = std::clock();
    pace::folded_matrix<uint32_t> folded_matrix(graph);
    end = std::clock();
    const double t_f = pace::test::time_in_ms(start, end);

    // test equality
    assert(pace::test::equals(ref_matrix, matrix, true));
    assert(pace::test::equals(ref_matrix, folded_matrix, true));

    const char fastest = t_f <= t_m ? 'f' : 'm';
    fmt::printf("%11s%11.3f%11.3f%11.3f%11c\n",  //
                filepath.filename(), t_r, t_m, t_f, fastest);
}

/**
 * @brief tests crossing_matrix with instances of dirpath
 */
void test_crossing_matrix_with(const fs::path dirpath) {
    fmt::printf("%11s%11s%11s%11s%11s\n",  //
                "instance", "reference", "matrix", "folded", "fastest");
    pace::test::print_line(67);
    for (const auto& file : std::filesystem::directory_iterator("medium_test_set/instances")) {
        if (!file.is_regular_file()) continue;
        test_crossing_matrix(file.path());
        std::cout << std::flush;
    }
}

/**
 * @brief tests matrix/matrix.hpp; its classes and functions
 */
int main(int argc, char** argv) {
    if (argc <= 1) {
        test_crossing_matrix_with("medium_test_set/instances");
    } else {
        const fs::path path = argv[1];
        if (path.extension() == ".gr") {
            test_crossing_matrix(path);
        } else {
            test_crossing_matrix(std::string{argv[1]} + "/instances");
        }
    }
    return 0;
}