#ifndef PACE2024_OUTPUT_HPP
#define PACE2024_OUTPUT_HPP

#include <vector>
#include <iostream>

namespace pace2024 {

template <typename T>
void print_output(const T n0, const std::vector<T>& ordering) {
    for (const T& element: ordering) {
        std::cout << element + n0 + 1 << std::endl;
    }
}

};  // namespace pace2024

#endif