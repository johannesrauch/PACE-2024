#ifndef PACE2024_TEST_UTILS_HPP
#define PACE2024_TEST_UTILS_HPP

#include <ctime>
#include <filesystem>
#include <iostream>

#include "printf.hpp"

namespace pace2024 {

namespace test {

/// @brief prints a line :o
void print_line(std::size_t length) {
    for (std::size_t i = 0; i < length; ++i) fmt::printf("-");
    fmt::printf("\n");
}

/// @brief returns the elapsed time between `start` and `end` in ms
double time_in_ms(const std::clock_t start, const std::clock_t end) {
    return 1000.0 * (end - start) / CLOCKS_PER_SEC;
}

template <typename R = uint32_t>
R get_ref_nof_crossings(std::filesystem::path filepath_instance) {
    auto filepath_nof_crossings = filepath_instance.parent_path().parent_path() /  //
                                  "nof_crossings" /                                //
                                  filepath_instance.filename();
    filepath_nof_crossings.replace_extension(".txt");
    // fmt::printf("%s\n", static_cast<std::string>(filepath_nof_crossings));
    std::ifstream file_nof_crossings(filepath_nof_crossings);
    assert(file_nof_crossings.good());
    R ref_nof_crossings;
    file_nof_crossings >> ref_nof_crossings;
    file_nof_crossings.close();
    return ref_nof_crossings;
}

};  // namespace test

};  // namespace pace2024

#endif