#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include "branch_and_cut.hpp"
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
    pace::branch_and_cut solver(instance);
    solver.run();

    uint32_t test = solver.get_nof_crossings();
    uint32_t ref = pace::test::get_ref_nof_crossings(instance.filepath);
    assert(test == ref);
    assert(test == pace::number_of_crossings(instance.graph(), solver.get_ordering()));
    (void) ref;

    fmt::printf("%s", instance.filepath);
}

void test_branch_and_cut(const fs::path dirpath) {
    for (const auto& file : fs::directory_iterator(dirpath)) {
        if (!file.is_regular_file()) continue;

        pace::instance instance(file.path());
        test_branch_and_cut_with(instance);
    }
}

/**
 * @brief tests branch_and_cut
 */
int main(int argc, char** argv) {
    if (argc <= 1) {
        test_branch_and_cut("tiny_test_set/instances");
        test_branch_and_cut("my_tests/instances");
    } else {
        test_branch_and_cut(std::string{argv[1]} + "/instances");
    }

    std::cout << "TEST::PACE::BRANCH_AND_CUT:\t\t\tOK" << std::endl;
    return 0;
}
