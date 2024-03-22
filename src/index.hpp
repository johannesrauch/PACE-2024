#ifndef PACE_INDEX_HPP
#define PACE_INDEX_HPP

#include <cassert>

namespace pace {

std::size_t flat_index(const std::size_t n, const std::size_t i, const std::size_t j) {
    assert(i < j);
    assert(i < n);
    assert(j < n);
    assert(n > 0);
    const std::size_t n_choose_2 = n * (n - 1) / 2;
    const std::size_t k = n_choose_2 - (n - i) * (n - i - 1) / 2 + j - i - 1;
    assert(k < n_choose_2);
    return k;
}

};  // namespace pace

#endif
