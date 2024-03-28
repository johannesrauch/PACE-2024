#include "branch_and_cut.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <string>

#include "crossings.hpp"
#include "debug_printf.hpp"
#include "input.hpp"
#include "instance.hpp"
#include "output.hpp"
#include "printf.hpp"
#include "test_utils.hpp"

namespace fs = std::filesystem;

template <typename T, typename R>
void test_branch_and_cut_with(const pace::instance<T, R>& instance) {
    const pace::bipartite_graph<T>& graph = instance.graph();
    fmt::printf("%11s%11u%11u%11u",            //
                instance.filepath.filename(),  //
                graph.get_n_fixed(),           //
                graph.get_n_free(),            //
                graph.get_m());
    std::cout << std::flush;

    pace::branch_and_cut solver(instance);
    std::clock_t start = std::clock();
    solver.template operator()<false>();
    std::clock_t end = std::clock();
    const double t = pace::test::time_in_ms(start, end);

    uint32_t test = solver.get_nof_crossings();
    (void)test;
    std::string warning;
    try {
        uint32_t ref = pace::test::get_ref_nof_crossings(instance.filepath);
        PACE_DEBUG_PRINTF("REF: %u\n\n", ref);
        assert(test == ref);
        assert(test == pace::number_of_crossings(instance.graph(), solver.get_ordering()));
        (void)ref;
    } catch (std::exception& e) {
        warning = e.what();
    }

    const pace::branch_and_cut_info& info = solver.get_info();
    fmt::printf("|%11.1f%11u%11u%11u|%11u%11u%11u|%11s\n",                             //
                t, info.nof_rows, info.nof_iterations, info.nof_branch_nodes,          //
                instance.get_lower_bound(), info.nof_crossings_h, info.nof_crossings,  //
                warning);
}

void test_branch_and_cut(const fs::path dirpath) {
    fmt::printf("%s\n\n", dirpath);
    fmt::printf("%11s%11s%11s%11s|%11s%11s%11s%11s|%11s%11s%11s|%11s\n",  //
                "instance", "n fixed", "n free", "m",                     //
                "time in ms", "nof rows", "nof iter", "nof nodes",        //
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
        pace::instance instance(filepath);
        test_branch_and_cut_with(instance);
    }
    fmt::printf("\n");
}

/**
 * @brief tests branch_and_cut
 */
int main(int argc, char** argv) {
    if (argc <= 1) {
        test_branch_and_cut("tiny_test_set/instances");
        test_branch_and_cut("my_tests/instances");
    } else {
        const fs::path path = argv[1];
        if (path.extension() == ".gr") {
            pace::instance instance(path);
            test_branch_and_cut_with(instance);
        } else {
            test_branch_and_cut(std::string{argv[1]} + "/instances");
        }
    }

    std::cout << "TEST::PACE::BRANCH_AND_CUT:\t\t\tOK" << std::endl;
    return 0;
}
