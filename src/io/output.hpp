#ifndef PACE_IO_OUTPUT_HPP
#define PACE_IO_OUTPUT_HPP

#include <vector>
#include <iostream>

namespace pace {

template <typename T>
void print_output(const std::size_t n_fixed, const std::vector<T>& ordering) {
    for (const T& element: ordering) {
        std::cout << element + n_fixed + 1 << std::endl;
    }
}

};  // namespace pace

#endif