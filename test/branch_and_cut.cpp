#include "branch_and_cut.hpp"
#include "output.hpp"

int main() {
    {
        pace2024::general_instance<uint16_t> instance(
            "tiny_test_set/matching_4_4.gr");
        pace2024::branch_and_cut<uint16_t> solver(instance);

        assert(solver.get_variable_index(0, 1) == 1);
        // assert(solver.get_variable_index(1, 0) == 1);
        assert(solver.get_variable_index(0, 3) == 3);
        // assert(solver.get_variable_index(3, 0) == 3);
        assert(solver.get_variable_index(1, 3) == 5);
        // assert(solver.get_variable_index(3, 1) == 5);
        assert(solver.get_variable_index(2, 3) == 6);
        // assert(solver.get_variable_index(3, 2) == 6);

        solver.solve(false);
    }

    pace2024::uint32_instance instance("my_tests/random_threshold.gr");
    pace2024::branch_and_cut<uint32_t> solver(instance);
    solver.solve();

    std::cout << "TEST::PACE2024::BRANCH_AND_CUT: OKAY" << std::endl;
    return 0;
}