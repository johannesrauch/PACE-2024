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

void test_highs_mip(pace::input& input) {
    pace::instance instance(input.get_graph());
    const pace::bipartite_graph& graph = instance.get_graph();
    fmt::printf("%11s%11u%11u%11u",         //
                input.filepath.filename(),  //
                graph.get_n_fixed(),        //
                graph.get_n_free(),         //
                graph.get_m());
    std::cout << std::flush;

    pace::branch_and_cut solver(instance);
    std::clock_t start = std::clock();
    solver();
    std::clock_t end = std::clock();
    const double t = pace::test::time_in_ms(start, end);

    crossing_number_t test = solver.upper_bound;
    (void)test;
    std::string warning;
    try {
        uint32_t ref = pace::test::get_ref_n_crossings(input.filepath);
        PACE_DEBUG_PRINTF("REF: %u\n\n", ref);
        assert(test == ref);
        assert(test == pace::number_of_crossings(instance.get_graph(), solver.get_ordering()));
        (void)ref;
    } catch (std::exception& e) {
        warning = e.what();
    }

    const pace::branch_and_cut_info& info = solver.get_info();
    fmt::printf("|%11.1f%11u%11u%11u|%11u%11u%11u|%11s\n",               //
                t, info.n_rows, info.n_iterations, info.n_branch_nodes,  //
                instance.get_lower_bound(), info.n_crossings_h, test,    //
                warning);
}

void test_highs_mip_with(const fs::path dirpath) {
    fmt::printf("%s\n\n", dirpath);
    fmt::printf("%11s%11s%11s%11s|%11s%11s%11s%11s|%11s%11s%11s|%11s\n",  //
                "instance", "n fixed", "n free", "m",                     //
                "time in ms", "n rows", "n iter", "n nodes",        //
                "lower bound", "heuristic", "optimal",                    //
                "warning");
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
 * @brief tests branch_and_cut
 */
int main(int argc, char** argv) {
    pace::test::pipe_cout_to_file raii_pipe("branch_and_cut.log");
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

    std::cout << "TEST::PACE::BRANCH_AND_CUT:\t\t\tOK" << std::endl;
    return 0;
}
