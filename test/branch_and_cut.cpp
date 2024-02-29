#include "branch_and_cut.hpp"

#include <filesystem>
#include <fstream>
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
    for (const auto& file : std::filesystem::directory_iterator("tiny_test_set")) {
        if (!file.is_regular_file()) continue;

        const auto filepath_instance = file.path();
        pace2024::uint32_bipartite_graph graph(filepath_instance);
        pace2024::branch_and_cut<uint32_t, uint32_t> solver(graph);
        solver.solve(false);

        uint32_t ref_nof_crossings;
        {
            auto filepath_nof_crossings =
                filepath_instance.parent_path() / "nof_crossings" / filepath_instance.filename();
            filepath_nof_crossings.replace_extension(".txt");
            // fmt::printf("%s\n", static_cast<std::string>(filepath_nof_crossings));
            std::ifstream file_nof_crossings(filepath_nof_crossings);
            assert(file_nof_crossings.good());
            file_nof_crossings >> ref_nof_crossings;
            file_nof_crossings.close();
        }

        assert(ref_nof_crossings == solver.get_nof_crossings());
        // fmt::printf("%s,%s\n", solver.get_nof_crossings(), solutions[i]);
        uint32_t test_nof_crossings =
            pace2024::crossing_number_of<uint32_t, uint32_t>(graph, solver.get_ordering());
        assert(ref_nof_crossings == test_nof_crossings);
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

    std::cout << "TEST::PACE2024::BRANCH_AND_CUT:\t\t\tOK" << std::endl;
    return 0;
}