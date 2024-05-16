#ifndef PACE_IO_OUTPUT_HPP
#define PACE_IO_OUTPUT_HPP

#include <iostream>
#include <vector>

namespace pace {

template <typename CharT, typename Traits, typename T, typename Allocator>
void print_output(  //
    std::basic_ostream<CharT, Traits>& ostream, const std::size_t n_fixed, const std::vector<T, Allocator>& ordering) {
    for (const T& element : ordering) {
        ostream << element + n_fixed + 1 << std::endl;
    }
}

template <typename T, typename Allocator>
void print_output(const std::size_t n_fixed, const std::vector<T, Allocator>& ordering) {
    print_output(std::cout, n_fixed, ordering);
}

};  // namespace pace

#endif