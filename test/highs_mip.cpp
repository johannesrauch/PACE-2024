#include "exact/highs_mip.hpp"

#include <set>

#include "exact/branch_and_cut.hpp"
#include "io/input.hpp"
#include "utils/test_utils.hpp"

namespace fs = std::filesystem;

using vertex_t = pace::vertex_t;
using crossing_number_t = pace::crossing_number_t;

void test_highs_mip(pace::input& input) {
    pace::instance instance(input.get_graph());
    const pace::bipartite_graph& graph = instance.get_graph();
    fmt::printf("%11s%11u%11u%11u",         //
                input.filepath.filename(),  //
                graph.get_n_fixed(),        //
                graph.get_n_free(),         //
                graph.get_m());
    std::cout << std::flush;

    pace::highs_mip solver(instance);
    std::vector<vertex_t> ordering;
    std::clock_t start = std::clock();
    solver(ordering);
    std::clock_t end = std::clock();
    const double t = pace::test::time_in_ms(start, end);

    crossing_number_t test = solver.upper_bound;
    PACE_DEBUG_PRINTF("TEST: %u\n", test);
    std::string warning;
    bool test_ok{true};
    try {
        crossing_number_t ref = pace::test::get_ref_n_crossings(input.filepath);
        PACE_DEBUG_PRINTF("REF:  %u\n\n", ref);
        test_ok &= test == ref;
        assert(test == ref);
        assert(test == pace::number_of_crossings(graph, ordering));
    } catch (std::exception& e) {
        warning = e.what();
    }

    fmt::printf("|%11.1f|%11u%11u|%11s%11s\n",             //
                t,                                         //
                solver.lower_bound(), solver.upper_bound,  //
                warning, test_ok ? "true" : "false");
}

void test_highs_mip_with(const fs::path dirpath) {
    fmt::printf("%s\n\n", dirpath);
    fmt::printf("%11s%11s%11s%11s|%11s|%11s%11s|%11s%11s\n",  //
                "instance", "n fixed", "n free", "m",         //
                "time in ms",                                 //
                "lb", "ub",                                   //
                "warning", "test ok");
    pace::test::print_line(136);
    std::set<fs::path> testcases;
    for (const auto& file : fs::directory_iterator(dirpath)) {
        if (file.is_regular_file()) {
            testcases.insert(file.path());
        }
    }
    for (const auto& filepath : testcases) {
        pace::input input(filepath);
        test_highs_mip(input);
    }
    fmt::printf("\n");
}

/**
 * @brief tests highs_mip
 */
int main(int argc, char** argv) {
    if (argc <= 1) {
        test_highs_mip_with("tiny_test_set/instances");
        test_highs_mip_with("my_tests/instances");
    } else {
        const fs::path path = argv[1];
        if (path.extension() == ".gr") {
            pace::input input(path);
            test_highs_mip(input);
        } else {
            test_highs_mip_with(std::string{argv[1]} + "/instances");
        }
    }
    return 0;
}
