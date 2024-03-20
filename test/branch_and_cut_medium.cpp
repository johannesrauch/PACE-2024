#include <filesystem>
#include <fstream>
#include <set>
#include <string>

#include "branch_and_cut.hpp"
#include "crossings.hpp"
#include "input.hpp"
#include "output.hpp"
#include "printf.hpp"
#include "test_utils.hpp"

void test_solver_w_instance(std::filesystem::path filepath_instance) {
    fmt::printf("%s\n", static_cast<std::string>(filepath_instance));
    std::cout << std::flush;

    pace::uint16_bipartite_graph graph;
    pace::parse_input(filepath_instance, graph);
    pace::branch_and_cut solver(graph);
    solver.solve(false);

    uint32_t ref_nof_crossings =
        pace::test::get_ref_nof_crossings<uint32_t>(filepath_instance);
    (void)ref_nof_crossings;

    uint32_t test_nof_crossings =
        pace::number_of_crossings<uint16_t, uint32_t>(graph, solver.get_ordering());
    (void)test_nof_crossings;

    PACE_DEBUG_PRINTF("%d=%d?\n", ref_nof_crossings, test_nof_crossings);
    assert(ref_nof_crossings == test_nof_crossings);
    assert(ref_nof_crossings == solver.get_nof_crossings());
}

/**
 * @brief test branch_and_cut solver on medium_test_set
 */
void test_solver_w_medium_test_set() {
    // sort first
    std::set<std::filesystem::path> sorted_directory;
    for (const auto& file : std::filesystem::directory_iterator("medium_test_set/instances")) {
        if (!file.is_regular_file()) continue;
        sorted_directory.insert(file.path());
    }

    // execute all testcases
    for (const auto& filepath_instance : sorted_directory) {
        test_solver_w_instance(filepath_instance);
    }
}

/**
 * @brief tests branch_and_cut on medium_test_set
 *
 * @return int
 */
int main(int argc, char** argv) {
    if (argc == 1) {
        test_solver_w_medium_test_set();
    } else {
        const std::string input = "medium_test_set/instances/" + std::string{argv[1]} + ".gr";
        test_solver_w_instance(input);
    }

    // test_solver_w_instance("medium_test_set/17.gr");
    std::cout << "TEST::PACE::BRANCH_AND_CUT_MEDIUM:\t\tOK" << std::endl;
    return 0;
}
