#include "branch_and_cut.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
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
    uint32_t ref = pace::test::get_ref_nof_crossings(instance.filepath);
    PACE_DEBUG_PRINTF("REFERENCE VALUE: %u\n", ref);
    assert(test == ref);
    assert(test == pace::number_of_crossings(instance.graph(), solver.get_ordering()));
    (void)ref;

    const pace::branch_and_cut_info& info = solver.get_info();
    fmt::printf("|%11.1f%11u%11u%11u\n", t, info.nof_rows, info.nof_iterations, info.nof_branch_nodes);
}

void test_branch_and_cut(const fs::path dirpath) {
    fmt::printf("%s\n\n", dirpath);
    fmt::printf("%11s%11s%11s%11s|%11s%11s%11s%11s\n",  //
                "instance", "n fixed", "n free", "m",   //
                "time in ms", "nof rows", "nof iter", "nof nodes");
    pace::test::print_line(90);
    for (const auto& file : fs::directory_iterator(dirpath)) {
        if (!file.is_regular_file()) continue;

        pace::instance instance(file.path());
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
