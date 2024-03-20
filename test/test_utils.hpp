#ifndef PACE2024_TEST_UTILS_HPP
#define PACE2024_TEST_UTILS_HPP

#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

#include "printf.hpp"

namespace pace2024 {

namespace test {

/// @brief prints a line :o
inline void print_line(std::size_t length) {
    for (std::size_t i = 0; i < length; ++i) fmt::printf("-");
    fmt::printf("\n");
}

/// @brief returns the elapsed time between `start` and `end` in ms
inline double time_in_ms(const std::clock_t start, const std::clock_t end) {
    return 1000.0 * (end - start) / CLOCKS_PER_SEC;
}

template <typename R = uint32_t>
R get_ref_nof_crossings(std::filesystem::path filepath_instance) {
    auto filepath_nof_crossings =                        //
        filepath_instance.parent_path().parent_path() /  //
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

template <typename T>
void get_ref_ordering(const std::filesystem::path filepath_instance,  //
                      const std::size_t n_fixed,                      //
                      const std::size_t n_free,                       //
                      std::vector<T> &ordering) {
    auto filepath_solution =                           //
        filepath_instance.parent_path().parent_path()  //
        / "solutions"                                  //
        / filepath_instance.filename();
    filepath_solution.replace_extension(".sol");
    std::ifstream file_solution(filepath_solution);
    assert(file_solution.good());

    ordering.clear();
    ordering.reserve(n_free);
    for (std::size_t i = 0; i < n_free; ++i) {
        T v;
        file_solution >> v;
        assert(v >= n_fixed + 1);
        v -= n_fixed + 1;
        assert(v < n_free);
        ordering.emplace_back(v);
    }
}

};  // namespace test

};  // namespace pace2024

#endif