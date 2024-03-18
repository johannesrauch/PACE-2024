#ifndef PACE2024_VECTOR_INTERSECTION_HPP
#define PACE2024_VECTOR_INTERSECTION_HPP

#include <cstddef>
#include <vector>

namespace pace2024 {

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
}

};  // namespace pace2024

#endif