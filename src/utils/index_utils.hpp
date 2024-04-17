#ifndef PACE_UTILS_INDEX_HPP
#define PACE_UTILS_INDEX_HPP

#include <cassert>
#include <tuple>

namespace pace {

/**
 * @brief maps the index pair (i, j), i < j, to a one-dimensional index
 *
 * @param n the minimum upper bound on a possible j
 * @param n_choose_2
 * @param i
 * @param j i < j
 * @return std::size_t one-dimensional index
 */
inline std::size_t  //
flat_index(const std::size_t n, const std::size_t n_choose_2, const std::size_t i, const std::size_t j) {
    assert(i < j);
    assert(i < n);
    assert(j < n);
    assert(n > 0);
    const std::size_t k = n_choose_2 - (n - i) * (n - i - 1) / 2 + j - i - 1;
    assert(k < n_choose_2);
    return k;
}

/**
 * @brief maps the index pair (i, j), i < j, to a one-dimensional index
 *
 * @param n the minimum upper bound on a possible j
 * @param i
 * @param j i < j
 * @return std::size_t one-dimensional index
 */
inline std::size_t flat_index(const std::size_t n, const std::size_t i, const std::size_t j) {
    return flat_index(n, n * (n - 1) / 2, i, j);
}

/**
 * @brief maps the pairs of indices of u < v < w to their one-dimensional indices
 *
 * @param n the minimum upper bound on a possible j
 * @param n_choose_2
 * @param u
 * @param v u < v
 * @param w v < w
 * @return std::tuple<std::size_t, std::size_t, std::size_t> indices of uv, vw, uw
 */
inline std::tuple<std::size_t, std::size_t, std::size_t>         //
flat_indices(const std::size_t n, const std::size_t n_choose_2,  //
             const std::size_t u, const std::size_t v, const std::size_t w) {
    const std::size_t uv = flat_index(n, n_choose_2, u, v);
    const std::size_t vw = flat_index(n, n_choose_2, v, w);
    const std::size_t uw = flat_index(n, n_choose_2, u, w);
    assert(uv < uw);
    assert(uw < vw);
    return std::make_tuple(uv, vw, uw);
}

};  // namespace pace

#endif
