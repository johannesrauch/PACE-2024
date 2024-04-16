#include "utils/crossings_utils.hpp"

#include "io/input.hpp"
#include "model/instance.hpp"
#include "utils/randomness_utils.hpp"
#include "utils/test_utils.hpp"

namespace fs = std::filesystem;

/**
 * @brief tests crossing_numbers_of function against crossing number matrix
 */
void test_crossing_numbers_of(const fs::path filepath) {
    pace::input input(filepath);
    pace::instance instance(input.get_graph());
    const pace::bipartite_graph& graph = instance.get_graph();
    const pace::crossing_matrix& cr_matrix = instance.get_cr_matrix();

    using vertex_t = pace::bipartite_graph::vertex_t;
    const std::size_t n_free = graph.get_n_free();
    for (vertex_t u = 0; u < n_free; ++u) {
        for (vertex_t v = u + 1; v < n_free; ++v) {
            auto [c_uv, c_vu] = pace::crossing_numbers_of(graph, u, v);
            assert(cr_matrix(u, v) == c_uv);
            assert(cr_matrix(v, u) == c_vu);
        }
    }
}

/**
 * @brief tests crossing_numbers_of with instances in dirpath
 */
void test_crossing_numbers_of_with(const fs::path dirpath) {
    fmt::printf("start test_crossing_numbers_of_with(%s)\n\n", dirpath);
    for (const auto& file : std::filesystem::directory_iterator(dirpath)) {
        if (!file.is_regular_file()) continue;
        test_crossing_numbers_of(file.path());
    }
    fmt::printf("end   test_crossing_numbers_of_with(%s)\n\n", dirpath);
}

/**
 * @brief tests the two number_of_crossings functions against each other
 */
void test_number_of_crossings(const fs::path filepath) {
    pace::input input(filepath);
    pace::instance instance(input.get_graph());
    const pace::bipartite_graph& graph = instance.get_graph();
    const pace::crossing_matrix& cr_matrix = instance.get_cr_matrix();

    using vertex_t = pace::bipartite_graph::vertex_t;
    std::vector<vertex_t> ordering(graph.get_n_free());
    pace::test::shuffle(ordering);

    std::clock_t start = std::clock();
    const uint32_t ref = pace::number_of_crossings(cr_matrix, ordering);
    std::clock_t end = std::clock();
    const double t_m = pace::test::time_in_ms(start, end);

    start = std::clock();
    const uint32_t test = pace::number_of_crossings(graph, ordering);
    end = std::clock();
    assert(ref == test);
    const double t_g = pace::test::time_in_ms(start, end);

    const char fastest = t_m <= t_g ? 'm' : 'g';
    fmt::printf("%11.3f%11.3f%11c\n", t_m, t_g, fastest);
}

/**
 * @brief tests number_of_crossings with instances in dirpath
 */
void test_number_of_crossings_with(const fs::path dirpath) {
    fmt::printf("start test_number_of_crossings_with(%s)\n\n", dirpath);
    fmt::printf("%11s%11s%11s%11s\n", "instance", "folded", "acctree", "fastest");
    pace::test::print_line(45);
    for (const auto& file : std::filesystem::directory_iterator(dirpath)) {
        if (!file.is_regular_file()) continue;
        test_number_of_crossings(file.path());
        std::cout << std::flush;
    }
    fmt::printf("end   test_number_of_crossings_with(%s)\n\n", dirpath);
}

/**
 * @brief tests utils/crossings_utils.hpp
 */
int main(int argc, char** argv) {
    if (argc == 1) {
        test_number_of_crossings_with("medium_test_set/instances");
        test_crossing_numbers_of_with("medium_test_set/instances");
    } else {
        const fs::path path = argv[1];
        if (path.extension() == ".gr") {
            test_number_of_crossings(path);
            test_crossing_numbers_of(path);
        } else {
            test_number_of_crossings_with(std::string{argv[1]} + "/instances");
            test_crossing_numbers_of_with(std::string{argv[1]} + "/instances");
        }
    }

    std::cout << "TEST::PACE::CROSSING_NUMBER:\t\tOK" << std::endl;
    return 0;
}