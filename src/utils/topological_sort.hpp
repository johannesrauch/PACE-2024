#ifndef PACE_UTILS_TOPOLOGICAL_SORT_HPP
#define PACE_UTILS_TOPOLOGICAL_SORT_HPP

#include <algorithm>
#include <cassert>
#include <stack>

#include "model/digraph.hpp"

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
 * @return true if a topological sort was found
 * @return false if a (directed) cycle was found
 */
template <bool RETURN_IF_CYCLIC = true, typename T, typename A>
bool topological_sort_dfs(const general_digraph<T>& graph,
                          std::vector<T, A>& ordering,       //
                          std::vector<uint8_t>& visited,  //
                          T v) {
    visited[v] = 1;
    auto& neighbors = graph.get_neighbors(v);
    bool acyclic = true;
    for (const T& u : neighbors) {
        if (visited[u] == 1 ||
            (visited[u] == 0 &&
             !topological_sort_dfs(graph, ordering, visited, u))) {
            // cycle found
            acyclic = false;
            if constexpr(RETURN_IF_CYCLIC) return acyclic;
        }
    }
    visited[v] = 2;
    ordering.emplace_back(v);
    return acyclic;
}

};

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
template <bool RETURN_IF_CYCLIC = true, typename T, typename A>
bool topological_sort(const general_digraph<T>& graph, std::vector<T, A>& ordering) {
    const std::size_t n = graph.get_n();
    ordering.clear();
    ordering.reserve(n);
    std::vector<uint8_t> visited(n, 0);

    bool acyclic = true;
    for (T v = 0; v < n; ++v) {
        if (visited[v] == 0 &&
            !internal::topological_sort_dfs<RETURN_IF_CYCLIC, T, A>(graph, ordering, visited, v)) {
            // cycle found
            acyclic = false;
            if constexpr(RETURN_IF_CYCLIC) return acyclic;
        }
    }

    std::reverse(ordering.begin(), ordering.end());
    assert(graph.get_n() == ordering.size());
    return acyclic;
}

};  // namespace pace

#endif