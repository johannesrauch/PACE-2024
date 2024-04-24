#include "exact/exact_solver.hpp"

#include <set>

#include "io/input.hpp"
#include "utils/test_utils.hpp"

namespace fs = std::filesystem;

using vertex_t = pace::vertex_t;
using crossing_number_t = pace::crossing_number_t;

void test_exact_solver(pace::input& input) {
    const pace::bipartite_graph& graph = input.get_graph();
    fmt::printf("%11s%11u%11u%11u",         //
                input.filepath.filename(),  //
                graph.get_n_fixed(),        //
                graph.get_n_free(),         //
                graph.get_m());
    std::cout << std::flush;

    pace::exact_solver solver(input);
    std::vector<vertex_t> ordering;
    std::clock_t start = std::clock();
    crossing_number_t test = solver(ordering);
    std::clock_t end = std::clock();
    const double t = pace::test::time_in_ms(start, end);

    PACE_DEBUG_PRINTF("TEST: %u\n", test);
    std::string warning;
    bool test_ok{true};
    try {
        uint32_t ref = pace::test::get_ref_n_crossings(input.filepath);
        PACE_DEBUG_PRINTF("REF:  %u\n\n", ref);
        test_ok &= test == ref;
        assert(test == ref);
        assert(test == pace::number_of_crossings(graph, ordering));
    } catch (std::exception& e) {
        warning = e.what();
    }

    fmt::printf("|%11.1f|%11s%11s\n", t, warning, test_ok ? "true" : "false");
}

void test_exact_solver_with(const fs::path dirpath) {
    fmt::printf("%s\n\n", dirpath);
    fmt::printf("%11s%11s%11s%11s|%11s|%11s%11s\n",    //
                "instance", "n fixed", "n free", "m",  //
                "time in ms", "warning", "test ok");
    pace::test::print_line(80);
    std::set<fs::path> testcases;
    for (const auto& file : fs::directory_iterator(dirpath)) {
        if (file.is_regular_file()) {
            testcases.insert(file.path());
        }
    }
    for (const auto& filepath : testcases) {
        pace::input input(filepath);
        test_exact_solver(input);
    }
    fmt::printf("\n");
}

/**
 * @brief tests branch_and_cut
 */
int main(int argc, char** argv) {
    if (argc <= 1) {
        test_exact_solver_with("tiny_test_set/instances");
        test_exact_solver_with("my_tests/instances");
    } else {
        const fs::path path = argv[1];
        if (path.extension() == ".gr") {
            pace::input input(path);
            test_exact_solver(input);
        } else {
            test_exact_solver_with(std::string{argv[1]} + "/instances");
        }
    }
    return 0;
}
