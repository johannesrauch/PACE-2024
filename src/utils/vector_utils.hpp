#ifndef PACE_UTILS_VECTOR_UTILS_HPP
#define PACE_UTILS_VECTOR_UTILS_HPP

#include <cstddef>
#include <vector>

#include "fmt/printf.hpp"

namespace pace {

//
// test functions
//

namespace test {

template <typename T, typename Allocator>
void print_vector(const std::vector<T, Allocator> &vec) {
    for (const T &elem : vec) {
        fmt::printf("%d, ", elem);
    }
    fmt::printf("\n");
}

template <typename T, class = typename std::enable_if_t<std::is_unsigned<T>::value>, typename Allocator>
bool is_permutation(const std::vector<T, Allocator> &vec) {
    const std::size_t n = vec.size();
    std::vector<bool> test(n);
    for (const T &v : vec) {
        if (v > n) return false;
        if (test[v]) return false;
        test[v] = true;
    }
    return true;
}

};  // namespace test

//
// vector utility functions
//

/**
 * @brief inverses a permutation of {0, 1, ..., in.size()}
 *
 * @param in permutation, in parameter
 * @param out out parameter
 */
template <typename T, typename Allocator1, typename Allocator2>
inline void inverse(const std::vector<T, Allocator1> &in, std::vector<T, Allocator2> &out) {
    assert(test::is_permutation(in));
    const std::size_t len = in.size();
    out.resize(len);
    for (std::size_t i = 0; i < len; ++i) {
        assert(0 <= in[i]);
        assert(in[i] < len);
        out[in[i]] = i;
    }
}

/**
 * @brief resizes and stores vec[0]=0, vec[1]=1, ..., vec[n-1]=n-1 in vec
 */
template <typename T, class = typename std::enable_if_t<std::is_unsigned<T>::value>, typename Allocator>
inline void identity(const std::size_t n, std::vector<T, Allocator> &vec) {
    vec.resize(n);
    for (std::size_t i = 0; i < n; ++i) vec[i] = i;
}

// todo: use stdlib functions

/**
 * @brief returns the median of the vector vec
 *
 * @tparam T element type of vec
 * @param vec vector
 * @return T the median of vec
 */
template <typename T, typename Allocator>
inline T median(const std::vector<T, Allocator> &vec) {
    const std::size_t len = vec.size();
    if (len == 0) {
        return 0;
    } else if (len == 1) {
        return vec[0];
    } else {
        const std::size_t len2 = len / 2;
        return len % 2 == 0 ? (vec[len2 - 1] + vec[len2]) / 2 : vec[len / 2];
    }
}

/**
 * @brief returns the number of common elements in the ascending ordered vectors
 *
 * @tparam T element type of vectors
 * @tparam R return type
 * @param vec1
 * @param vec2
 * @return R number of common elements
 */
template <typename T, typename R = std::size_t>
R sorted_vector_intersection(const std::vector<T> &vec1, const std::vector<T> &vec2) {
    R size_of_intersection{0};
    for (std::size_t i = 0, j = 0; i < vec1.size() && j < vec2.size();) {
        if (vec1[i] == vec2[j]) {
            ++size_of_intersection;
            ++i;
            ++j;
        } else if (vec1[i] < vec2[j]) {
            ++i;
        } else {
            ++j;
        }
    }
    return size_of_intersection;
}

/**
 * @brief returns the union of the ascending sorted vectors `vec1` and `vec2` in `out`
 *
 * @tparam T element type of vectors
 * @param vec1
 * @param vec2
 * @param out out-parameter, union of `vec1` and `vec2` is stored here
 */
template <typename T>
void sorted_vector_union(const std::vector<T> &vec1, const std::vector<T> &vec2, std::vector<T> &out) {
    out.clear();
    out.reserve(vec1.size() + vec2.size());
    std::size_t i = 0, j = 0;
    std::size_t len1 = vec1.size(), len2 = vec2.size();
    while (i < len1 && j < len2) {
        if (vec1[i] == vec2[j]) {
            out.emplace_back(vec1[i]);
            ++i;
            ++j;
        } else if (vec1[i] < vec2[j]) {
            out.emplace_back(vec1[i]);
            ++i;
        } else {
            out.emplace_back(vec2[j]);
            ++j;
        }
    }
    while (i < len1) out.emplace_back(vec1[i++]);
    while (j < len2) out.emplace_back(vec2[j++]);
}

};  // namespace pace

#endif