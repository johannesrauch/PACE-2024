#ifndef PACE2024_TOPOLOGICAL_SORT_HPP
#define PACE2024_TOPOLOGICAL_SORT_HPP

#include <algorithm>
#include <cassert>
#include <stack>

#include "graph.hpp"

namespace pace2024 {

/**
 * @brief a depth-first search through the given graph.
 * used by topological_sort below.
 *
 * @tparam T vertex type
 * @param graph
 * @param ordering vector, where is stored
 * @param visited
 * @param v current vertex
 * @return true if a topological sort was found
 * @return false if a (directed) cycle was found
 */
template <typename T>
bool topological_sort_dfs(const general_digraph<T>& graph,
                          std::vector<T>& ordering,       //
                          std::vector<uint8_t>& visited,  //
                          T v) {
    visited[v] = 1;
    auto& adjacency_list = graph.get_adjacency_list(v);
    for (T u : adjacency_list) {
        if (visited[u] == 1 ||
            (visited[u] == 0 &&
             !topological_sort_dfs(graph, ordering, visited, u))) {
            // cycle found
            return false;
        }
    }
    visited[v] = 2;
    ordering.emplace_back(v);
    return true;
}

/**
 * @brief computes a topological sort (if possible) of the given graph
 * and stores it in the vector ordering.
 *
 * @tparam T vertex type
 * @param graph
 * @param ordering
 * @return true if a topological sort was found
 * @return false if a (directed) cycle was found
 */
template <typename T>
bool topological_sort(const general_digraph<T>& graph, std::vector<T>& ordering) {
    const T n = graph.get_n();
    ordering.clear();
    std::vector<uint8_t> visited(n, 0);

    for (T v = 0; v < n; ++v) {
        if (visited[v] == 0 &&
            !topological_sort_dfs(graph, ordering, visited, v)) {
            // cycle found
            return false;
        }
    }

    std::reverse(ordering.begin(), ordering.end());
    assert(graph.get_n() == ordering.size());
    return true;
}

};  // namespace pace2024

#endif