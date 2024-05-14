#include "exact/branch_and_cut.hpp"

#include <set>

#include "io/input.hpp"
#include "io/output.hpp"
#include "model/instance.hpp"
#include "utils/crossings_utils.hpp"
#include "utils/test_utils.hpp"

namespace fs = std::filesystem;

using vertex_t = pace::vertex_t;
using crossing_number_t = pace::crossing_number_t;

void test_branch_and_cut(pace::input& input) {
    pace::instance instance(input.get_graph());
    const pace::bipartite_graph& graph = instance.get_graph();
    fmt::printf("%11s%11u%11u%11u",         //
                input.filepath.filename(),  //
                graph.get_n_fixed(),        //
                graph.get_n_free(),         //
                graph.get_m());
    std::cout << std::flush;

    pace::branch_and_cut solver(instance);
    std::vector<vertex_t> ordering;
    const std::chrono::time_point<std::chrono::system_clock> t0 = std::chrono::system_clock::now();
    crossing_number_t test = solver(ordering);

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

    const pace::branch_and_cut_info& info = solver.get_info();
    fmt::printf("|%11.1f%11u%11u%11u|%11u%11u%11u|%11s%11s\n",                                         //
                pace::elapsed_walltime_in_s(t0), info.n_rows, info.n_iterations, info.n_branch_nodes,  //
                info.lower_bound, info.n_crossings_h, info.upper_bound,                                //
                warning, test_ok ? "true" : "false");
}

void test_branch_and_cut_with(const fs::path dirpath) {
    fmt::printf("%s\n\n", dirpath);
    fmt::printf("%11s%11s%11s%11s|%11s%11s%11s%11s|%11s%11s%11s|%11s%11s\n",  //
                "instance", "n fixed", "n free", "m",                         //
                "time in s", "n rows", "n iter", "n nodes",                   //
                "lower bound", "heuristic", "optimal",                        //
                "warning", "test ok");
    pace::test::print_line(147);
    std::set<fs::path> testcases;
    for (const auto& file : fs::directory_iterator(dirpath)) {
        if (file.is_regular_file()) {
            testcases.insert(file.path());
        }
    }
    for (const auto& filepath : testcases) {
        pace::input input(filepath);
        test_branch_and_cut(input);
    }
    fmt::printf("\n");
}

/**
 * @brief tests branch_and_cut
 */
int main(int argc, char** argv) {
    if (argc <= 1) {
        test_branch_and_cut_with("tiny_test_set/instances");
        test_branch_and_cut_with("my_tests/instances");
    } else {
        const fs::path path = argv[1];
        if (path.extension() == ".gr") {
            pace::input input(path);
            test_branch_and_cut(input);
        } else {
            test_branch_and_cut_with(std::string{argv[1]} + "/instances");
        }
    }
    return 0;
}
