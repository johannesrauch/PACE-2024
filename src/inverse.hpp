#ifndef PACE2024_INVERSE_HPP
#define PACE2024_INVERSE_HPP

#include <vector>

namespace pace2024 {

/**
 * @brief inverses a permutation of {0, 1, ..., in.size()}
 * 
 * @param in permutation, in parameter
 * @param out out parameter
 */
template <typename T>
inline void inverse(const std::vector<T> &in, std::vector<T> &out) {
    const std::size_t len = in.size();
    out.resize(len);
    for (std::size_t i = 0; i < len; ++i) {
        assert(0 <= in[i]);
        assert(in[i] < len);
        out[in[i]] = i;
    }
}

};  // namespace pace2024

#endif
