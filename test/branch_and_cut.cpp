#include "branch_and_cut.hpp"

#include <filesystem>
#include <fstream>
#include <string>

#include "input.hpp"
#include "crossing_number.hpp"
#include "output.hpp"
#include "printf.hpp"

namespace fs = std::filesystem;

/**
 * @brief test branch_and_cut solver on tiny_test_set
 */
void test_solver_w_tiny_test_set() {
    for (const auto& file : fs::directory_iterator("tiny_test_set")) {
        if (!file.is_regular_file()) continue;

        const auto filepath_instance = file.path();
        pace2024::uint32_bipartite_graph graph;
        pace2024::parse_input(filepath_instance, graph);
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
            pace2024::number_of_crossings<uint32_t, uint32_t>(graph, solver.get_ordering());
        assert(ref_nof_crossings == test_nof_crossings);
    }
}

/**
 * @brief test branch_and_cut solver on a random instance
 */
void test_solver_w_random_instance() {
    const fs::path filepath = "my_tests/random_threshold.gr";
    pace2024::uint32_bipartite_graph graph;
    pace2024::parse_input(filepath, graph);
    pace2024::branch_and_cut<uint32_t, uint32_t> solver(graph);
    solver.solve(false);
    uint32_t nof_crossings = solver.get_nof_crossings();
    // std::cout << nof_crossings << std::endl;
    assert(nof_crossings == 414);
    uint32_t nof_crossings_2 =
        pace2024::number_of_crossings<uint32_t, uint32_t>(graph, solver.get_ordering());
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
