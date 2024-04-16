#ifndef PACE_OUTPUT_HPP
#define PACE_OUTPUT_HPP

#include <vector>
#include <iostream>

namespace pace {

template <typename T>
void print_output(const std::size_t n0, const std::vector<T>& ordering) {
    for (const T& element: ordering) {
        std::cout << element + n0 + 1 << std::endl;
    }
}

};  // namespace pace

#endif