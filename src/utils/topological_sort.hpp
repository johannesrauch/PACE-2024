#ifndef PACE_UTILS_TOPOLOGICAL_SORT_HPP
#define PACE_UTILS_TOPOLOGICAL_SORT_HPP

#include <algorithm>
#include <cassert>
#include <stack>

#include "model/digraph.hpp"
#include "utils/randomness_utils.hpp"

namespace pace {

namespace internal {

/**
 * @brief a depth-first search through the given graph.
 * used by topological_sort below.
 *
 * @tparam T vertex type
 * @tparam A allocator
 * @param graph
 * @param ordering vector, where is stored
 * @param visited
 * @param v current vertex
 * @return number of cycles found
 */
template <bool RETURN_IF_CYCLIC = true, typename T, typename A>
std::size_t topological_sort_dfs(const general_digraph<T>& graph,  //
                                 std::vector<T, A>& ordering,      //
                                 std::vector<uint8_t>& visited,    //
                                 T v, std::vector<std::pair<T, T>>& backarcs) {
    visited[v] = 1;
    const auto& neighbors = graph.get_neighbors(v);
    std::size_t n_cycles = 0;
    for (const T& u : neighbors) {
        std::size_t n_cycles_ = 0;
        if (visited[u] == 0) {
            n_cycles_ = topological_sort_dfs<RETURN_IF_CYCLIC, T, A>(
                graph, ordering, visited, u, backarcs);
        } else if (visited[u] == 1) {
            n_cycles_ = 1;
            backarcs.emplace_back(v, u);
        }
        n_cycles += n_cycles_;
        if (RETURN_IF_CYCLIC && n_cycles > 0) {
            return 1;
        }
    }
    visited[v] = 2;
    ordering.emplace_back(v);
    return n_cycles;
}

};  // namespace internal

/**
 * @brief computes a topological sort (if possible) of the given graph
 * and stores it in the vector ordering.
 *
 * @tparam T vertex type
 * @tparam A allocator
 * @param graph
 * @param ordering
 * @return true if a topological sort was found
 * @return false if a (directed) cycle was found
 */
template <typename T, typename A>
bool topological_sort(const general_digraph<T>& graph,
                      std::vector<T, A>& ordering) {
    const std::size_t n = graph.get_n();
    ordering.clear();
    ordering.reserve(n);
    std::vector<uint8_t> visited(n, 0);
    std::vector<std::pair<T, T>> backarcs;

    for (T v = 0; v < n; ++v) {
        if (visited[v] == 0 && internal::topological_sort_dfs(
                                   graph, ordering, visited, v, backarcs) > 0) {
            // cycle found
            return false;
        }
    }

    std::reverse(ordering.begin(), ordering.end());
    assert(n == ordering.size());
    return true;
}

/**
 * @brief computes a topological sort (if possible) of the given graph
 * and stores it in the vector ordering.
 *
 * @tparam T vertex type
 * @tparam A allocator
 * @param graph
 * @param ordering
 * @return number of cycles found
 */
template <typename T, typename A>
std::size_t topological_sort_rd(const general_digraph<T>& graph,
                                std::vector<T, A>& ordering) {
    const std::size_t n = graph.get_n();
    ordering.clear();
    ordering.reserve(n);
    std::vector<uint8_t> visited(n, 0);
    std::vector<T> vertices(n);
    pace::test::shuffle(vertices);  // may useful for round heuristic
    std::vector<std::pair<T, T>> backarcs;

    std::size_t n_cycles = 0;
    for (const T& v : vertices) {
        if (visited[v] == 0) {
            n_cycles += internal::topological_sort_dfs<false, T, A>(
                graph, ordering, visited, v, backarcs);
        }
    }

    std::reverse(ordering.begin(), ordering.end());
    assert(n == ordering.size());
    return n_cycles;
}

/**
 * @brief computes a topological sort (if possible) of the given graph
 * and stores it in the vector ordering.
 *
 * @tparam T vertex type
 * @tparam A allocator
 * @param graph
 * @param ordering
 * @param backarcs
 * @return number of cycles found
 */
template <typename T, typename A>
std::size_t topological_sort_backarcs(const general_digraph<T>& graph,
                                      std::vector<T, A>& ordering,
                                      std::vector<std::pair<T, T>>& backarcs) {
    const std::size_t n = graph.get_n();
    ordering.clear();
    ordering.reserve(n);
    std::vector<uint8_t> visited(n, 0);

    std::size_t n_cycles = 0;
    for (T v = 0; v < n; ++v) {
        if (visited[v] == 0) {
            n_cycles += internal::topological_sort_dfs<false, T, A>(
                graph, ordering, visited, v, backarcs);
        }
    }

    std::reverse(ordering.begin(), ordering.end());
    assert(n == ordering.size());
    return n_cycles;
}

};  // namespace pace

#endif