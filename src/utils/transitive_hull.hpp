#ifndef PACE_UTILS_TRANSITIVE_HULL_HPP
#define PACE_UTILS_TRANSITIVE_HULL_HPP

#include <unordered_set>

#include "model/digraph.hpp"
#include "utils/index_utils.hpp"
#include "utils/topological_sort.hpp"

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
void transitive_hull_dfs(const general_digraph<T> &graph, std::vector<uint8_t> &visited, const T u) {
    visited[u] += 2;
    for (const T &v : graph.get_neighbors(u)) {
        if (visited[v] <= 1) {
            transitive_hull_dfs(graph, visited, v);
        }
    }
}

};  // namespace internal

/**
 * @brief computes the transitive hull of digraph
 *
 * @tparam T vertex type
 * @param graph input digraph
 * @param new_arcs out parameter; new arcs are stored here
 */
template <typename T>
void transitive_hull(const general_digraph<T> &graph, std::vector<std::pair<T, T>> &new_arcs) {
    const std::size_t n = graph.get_n();
    std::vector<uint8_t> visited(n);

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
}

// todo: speed the acyclic method up

/**
 * @brief computes the transitive hull of an cyclic digraph (Purdom's algorithm)
 *
 * @tparam T vertex type
 * @param graph input digraph
 * @param new_arcs out parameter; new arcs are stored here
 */
template <typename T>
void transitive_hull_of_acyclic(const general_digraph<T> &graph, std::vector<std::pair<T, T>> &new_arcs) {
    const std::size_t n = graph.get_n();
    assert(n > 0);
    std::vector<T> ordering;
    bool is_acyclic = topological_sort(graph, ordering);
    (void)is_acyclic;  // to suppress warning
    assert(is_acyclic);

    std::vector<std::unordered_set<T>> adjacency_lists(n);
    for (std::size_t u = 0; u < n; ++u) {
        const std::vector<T> &nbors = graph.get_neighbors(u);
        adjacency_lists[u].insert(nbors.begin(), nbors.end());
    }

    new_arcs.clear();
    for (std::size_t i = n - 1;; --i) {
        const T &u = ordering[i];
        std::unordered_set<T> new_neighbors;
        for (const T &v : adjacency_lists[u]) {
            new_neighbors.insert(adjacency_lists[v].begin(), adjacency_lists[v].end());
        }
        for (const T &v : new_neighbors) {
            if (adjacency_lists[u].count(v) == 0) {
                new_arcs.emplace_back(u, v);
            }
        }
        adjacency_lists[u].insert(new_neighbors.begin(), new_neighbors.end());

        if (i == 0) break;
    }
}

};  // namespace pace

#endif
