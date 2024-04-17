#include "median_heuristic.hpp"

#include <cassert>
#include <filesystem>
#include <iostream>

#include "crossings.hpp"
#include "instance.hpp"
#include "printf.hpp"

/**
 * @brief tests the heuristic solvers median_heuristic and probmedian_heuristic
 * (and shift_heuristic in the background)
 * 
 * goals:
 * - no compilation error
 * - no crash
 * - probmedian should be better
 */
template <typename T, typename R>
void test_median_heuristics(const pace::instance<T, R>& instance) {
    std::vector<T> ordering;
    uint32_t n_crossings = pace::median_heuristic{instance}(ordering);
    uint32_t nof_crossings_ = pace::probmedian_heuristic{instance}(ordering);
    assert(nof_crossings_ <= n_crossings);
}

/**
 * @brief tests median_heuristic.hpp
 */
int main() {
    for (const auto& file : std::filesystem::directory_iterator("tiny_test_set/instances")) {
        if (!file.is_regular_file()) continue;
        fmt::printf("%s\n", file.path());

        pace::instance instance(file.path());
        test_median_heuristics(instance);
    }

    std::cout << "TEST::PACE::MEDIAN_HEURISTIC:\t\tOK" << std::endl;
    return 0;
}