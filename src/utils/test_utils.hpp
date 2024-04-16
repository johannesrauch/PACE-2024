#ifndef PACE_UTILS_TEST_UTILS_HPP
#define PACE_UTILS_TEST_UTILS_HPP

#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

#include "fmt/printf.hpp"

namespace pace {

namespace test {

/**
 * @brief prints a line :o
 */
void print_line(std::size_t length) {
    for (std::size_t i = 0; i < length; ++i) fmt::printf("-");
    fmt::printf("\n");
}

/**
 * @brief returns the elapsed time between start and end in ms
 */
inline double time_in_ms(const std::clock_t start, const std::clock_t end) {
    return 1000.0 * (end - start) / CLOCKS_PER_SEC;
}

/**
 * @brief loads and returns reference number of crossings
 * 
 * @param filepath_instance filepath to instance (!)
 */
uint32_t get_ref_n_crossings(std::filesystem::path filepath_instance) {
    auto filepath_n_crossings =                        //
        filepath_instance.parent_path().parent_path() /  //
        "nof_crossings" /                                //
        filepath_instance.filename();
    filepath_n_crossings.replace_extension(".txt");

    std::ifstream file_n_crossings(filepath_n_crossings);
    if (!file_n_crossings.good()) throw std::runtime_error("no ref");
    
    uint32_t ref_n_crossings;
    file_n_crossings >> ref_n_crossings;
    file_n_crossings.close();
    return ref_n_crossings;
}

/**
 * @brief loads reference ordering from filepath
 * 
 * @tparam T vertex type
 * @param filepath_instance filepath to instance (!)
 */
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
    if (!file_solution.good()) throw std::runtime_error("no ref");

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

};  // namespace pace

#endif