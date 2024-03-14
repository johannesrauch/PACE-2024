#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include "branch_and_cut.hpp"
#include "crossings.hpp"
#include "debug_printf.hpp"
#include "input.hpp"
#include "output.hpp"
#include "printf.hpp"
#include "test_utils.hpp"

namespace fs = std::filesystem;

/**
 * @brief test branch_and_cut solver on tiny_test_set
 */
void test_solver_w_tiny_test_set() {
    for (const auto& file : fs::directory_iterator("tiny_test_set/instances")) {
        if (!file.is_regular_file()) continue;

        const auto filepath_instance = file.path();
        std::cout << filepath_instance << std::endl
                  << std::flush;
        pace2024::bipartite_graph graph;
        pace2024::parse_input(filepath_instance, graph);
        pace2024::branch_and_cut solver(graph);
        solver.solve(false);

        uint32_t ref_nof_crossings = pace2024::test::get_ref_nof_crossings(filepath_instance);
        PACE2024_DEBUG_PRINTF("%u=%u?\n", ref_nof_crossings, solver.get_nof_crossings());
        assert(ref_nof_crossings == solver.get_nof_crossings());
        // fmt::printf("%s,%s\n", solver.get_nof_crossings(), solutions[i]);
        uint32_t test_nof_crossings = pace2024::number_of_crossings(graph, solver.get_ordering());
        assert(ref_nof_crossings == test_nof_crossings);
    }
}

/**
 * @brief test branch_and_cut solver on a random instance
 */
void test_solver_w_random_instance() {
    const fs::path filepath = "my_tests/random_threshold.gr";
    std::cout << filepath << std::endl
              << std::flush;
    pace2024::bipartite_graph graph;
    pace2024::parse_input(filepath, graph);
    pace2024::branch_and_cut solver(graph);
    solver.solve(false);
    uint32_t nof_crossings = solver.get_nof_crossings();
    // std::cout << nof_crossings << std::endl;
    assert(nof_crossings == 414);
    uint32_t nof_crossings_2 = pace2024::number_of_crossings(graph, solver.get_ordering());
    assert(nof_crossings == nof_crossings_2);
}

/**
 * @brief tests branch_and_cut
 */
int main() {
    test_solver_w_tiny_test_set();
    test_solver_w_random_instance();

    std::cout << "TEST::PACE2024::BRANCH_AND_CUT:\t\t\tOK" << std::endl;
    return 0;
}
