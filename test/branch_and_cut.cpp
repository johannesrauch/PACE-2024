#include "branch_and_cut.hpp"

#include <string>

#include "crossing_number.hpp"
#include "output.hpp"
#include "printf.hpp"

/**
 * @brief tests get_variable_index of branch_and_cut
 */
void test_get_variable_index() {
    const std::string filepath = "tiny_test_set/matching_4_4.gr";
    pace2024::uint32_bipartite_graph graph(filepath);
    pace2024::branch_and_cut<uint32_t, uint16_t> solver(graph);

    assert(solver.get_variable_index(0, 1) == 1);
    // assert(solver.get_variable_index(1, 0) == 1);
    assert(solver.get_variable_index(0, 3) == 3);
    // assert(solver.get_variable_index(3, 0) == 3);
    assert(solver.get_variable_index(1, 3) == 5);
    // assert(solver.get_variable_index(3, 1) == 5);
    assert(solver.get_variable_index(2, 3) == 6);
    // assert(solver.get_variable_index(3, 2) == 6);
}

/**
 * @brief test branch_and_cut solver on tiny_test_set
 */
void test_solver_w_tiny_test_set() {
    const std::string tiny_tests[13] =
        {"complete_4_5.gr",
         "cycle_8_shuffled.gr",
         "cycle_8_sorted.gr",
         "grid_9_shuffled.gr",
         "ladder_4_4_shuffled.gr",
         "ladder_4_4_sorted.gr",
         "matching_4_4.gr",
         "path_9_shuffled.gr",
         "path_9_sorted.gr",
         "plane_5_6.gr",
         "star_6.gr",
         "tree_6_10.gr",
         "website_20.gr"};
    const uint16_t solutions[13] = {60, 4, 3, 17, 11, 3, 0, 6, 0, 0, 0, 13, 17};

    for (uint16_t i = 0; i < 13; ++i) {
        const std::string filepath = "tiny_test_set/" + tiny_tests[i];
        pace2024::uint32_bipartite_graph graph(filepath);
        pace2024::branch_and_cut<uint32_t, uint32_t> solver(graph);
        solver.solve(false);
        // assert(solver.get_nof_crossings() == solutions[i]);
        fmt::printf("%s,%s\n", solver.get_nof_crossings(), solutions[i]);
        std::cout.flush();
    }
}

/**
 * @brief test branch_and_cut solver on a random instance
 */
void test_solver_w_random_instance() {
    const std::string filepath = "my_tests/random_threshold.gr";
    pace2024::uint32_bipartite_graph graph(filepath);
    pace2024::branch_and_cut<uint32_t, uint32_t> solver(graph);
    solver.solve(false);
    uint32_t nof_crossings = solver.get_nof_crossings();
    // std::cout << nof_crossings << std::endl;
    assert(nof_crossings == 414);
    uint32_t nof_crossings_2 =
        pace2024::crossing_number_of<uint32_t, uint32_t>(graph, solver.get_ordering());
    assert(nof_crossings == nof_crossings_2);
}

/**
 * @brief tests branch_and_cut
 * 
 * @return int 
 */
int main() {
    test_get_variable_index();
    test_solver_w_tiny_test_set();
    test_solver_w_random_instance();

    std::cout << "TEST::PACE2024::BRANCH_AND_CUT: OKAY" << std::endl;
    return 0;
}