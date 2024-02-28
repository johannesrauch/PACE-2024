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
R vector_intersection(const std::vector<T>& vec1, const std::vector<T>& vec2) {
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

};

#endif