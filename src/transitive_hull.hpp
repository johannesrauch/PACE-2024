#ifndef PACE_TRANSITIVE_HULL_HPP
#define PACE_TRANSITIVE_HULL_HPP

#include "digraph.hpp"

namespace pace {

namespace internal {

/**
 * @brief internal dfs method for transitive_hull
 * 
 * @tparam T vertex type
 * @param graph input digraph
 * @param visited visited[u] = 0 (not visited), = 1 (neighbor of s), >= 2 (visited)
 * @param u current vertex
 */
template <typename T>
void transitive_hull_dfs(const digraph<T> &graph, std::vector<uint8_t> &visited, const T u) {
    visited[v] += 2;
    for (const T& v : graph.get_neighbors(u)) {
        if (visited[v] <= 1) {
            transitive_hull_dfs(graph, visited, v);
        }
    }
}

};  // namespace internal

/**
 * @brief computes the transitive hull of graph
 * 
 * @tparam T vertex type
 * @param graph input digraph
 */
template <typename T>
void transitive_hull(digraph<T> &graph) {
    const std::size_t n = graph.get_n();
    std::vector<uint8_t> visited(n);
    std::vector<std::pair<T, T>> new_arcs;

    for (T s = 0; s < n; ++s) {
        // initialize visited, visited[v] = 1 if v is a nbor of s or s and 0 otherwise
        for (uint8_t &vis : visited) vis = 0;
        visited[s] = 1;
        for (const T &v : graph.get_neighbors(s)) visited[v] = 1;

        internal::transitive_hull_dfs(graph, visited, s);

        // store new edges (for performance we don't include them right away)
        for (T v = 0; v < n; ++v) {
            if (visited[v] == 2) {
                new_arcs.emplace_back(s, v);
            }
        }
    }

    for (const auto &[u, v] : new_arcs) {
        graph.add_arc(u, v);
    }
}

};  // namespace pace

#endif
