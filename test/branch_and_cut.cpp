#include "branch_and_cut.hpp"

int main() {
    pace2024::instance instance("tiny_test_set/matching_4_4.gr");
    pace2024::branch_and_cut<uint16_t> solver(instance);

    assert(solver.get_variable_index(0, 1) == 1);
    assert(solver.get_variable_index(1, 0) == 1);
    assert(solver.get_variable_index(0, 3) == 3);
    assert(solver.get_variable_index(3, 0) == 3);
    assert(solver.get_variable_index(1, 3) == 5);
    assert(solver.get_variable_index(3, 1) == 5);
    assert(solver.get_variable_index(2, 3) == 6);
    assert(solver.get_variable_index(3, 2) == 6);

    std::cout << "TEST::PACE2024::BRANCH_AND_CUT: OKAY" << std::endl;
    return 0;
}