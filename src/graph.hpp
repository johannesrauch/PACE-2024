#ifndef PACE2024_GRAPH_HPP
#define PACE2024_GRAPH_HPP

#include <cstdint>
#include <queue>
#include <vector>

#include "debug_printf.hpp"

namespace pace2024 {

/**
 * @brief class for a generic directed graph with
 * definable vertex type
 *
 * @tparam T vertex type
 */
template <typename T>
class general_digraph {
   private:
    /// @brief stores all neighbors of each vertex
    std::vector<std::vector<T>> adjacency_lists;

   public:
    using vertextype = T;

    // delete copy constructor, move constructor, copy assignment and move assignment
    general_digraph(const general_digraph &rhs) = delete;
    general_digraph(general_digraph &&rhs) = delete;
    general_digraph &operator=(const general_digraph &rhs) = delete;
    general_digraph &operator=(general_digraph &&rhs) = delete;

    /**
     * @brief constructs an empty digraph with n vertices
     *
     * @param n number of vertices
     */
    general_digraph(const std::size_t n)
        : adjacency_lists(n),
          nof_neighbors_in_rollback_base(n) {}

    /**
     * @brief add arc (u, v)
     *
     * @param u vertex, 0 <= u < get_n()
     * @param v vertex, 0 <= v < get_n()
     */
    void add_arc(const T u, const T v) {
        assert(u < get_n());
        assert(v < get_n());
        adjacency_lists[u].emplace_back(v);
    }

    /**
     * @brief return number of vertices
     *
     * @return std::size_t number of vertices
     */
    std::size_t get_n() const {
        return adjacency_lists.size();
    }

    /**
     * @brief get neighbors of vertex v
     *
     * @param v vertex, 0 <= v < get_n()
     * @return const std::vector<T>& neighbors of v
     */
    const std::vector<T> &get_adjacency_list(const T v) const {
        return adjacency_lists[v];
    }

    /// @brief deletes all arcs
    void delete_arcs() {
        for (std::size_t i = 0; i < get_n(); ++i) {
            adjacency_lists[i].clear();
        }
    }

    /**
     * @brief checks if u and v are connected by a directed path
     * with a breadth-first search
     *
     * @param u vertex, 0 <= u < get_n()
     * @param v vertex, 0 <= v < get_n()
     * @return true if u and v are connected by a directed path
     * @return false otherwise
     */
    bool connected(const T u, const T v) {
        std::queue<T> q;
        q.emplace(u);
        std::vector<bool> visited(get_n(), false);

        while (!q.empty()) {
            const T current = q.front();
            q.pop();

            if (current == v) return true;
            for (const T &neighbor : adjacency_lists[current]) {
                if (!visited[neighbor]) q.emplace(neighbor);
            }
        }
        return false;
    }
};

namespace test {

/**
 * @brief tests if the given graph is a tournament
 * (a tournament is a directed graph such that (i, j) or (j, i) is an arc for every i < j)
 * not efficient!
 *
 * @tparam T vertex type
 * @param graph
 * @return true if the graph is a tournament
 * @return false otherwise
 */
template <typename T>
bool is_tournament(const general_digraph<T> &graph) {
    const std::size_t n = graph.get_n();
    for (T i = 0; i < n - 1; ++i) {
        const auto &neighbors_i = graph.get_adjacency_list(i);
        for (T j = i + 1; j < n; ++j) {
            const auto &neighbors_j = graph.get_adjacency_list(j);
            if (std::find(neighbors_i.begin(), neighbors_i.end(), j) == neighbors_i.end() &&
                std::find(neighbors_j.begin(), neighbors_j.end(), i) == neighbors_j.end()) {
                PACE2024_DEBUG_PRINTF("%s,%s\n", i, j);
                return false;
            }
        }
    }
    return true;
}

};  // namespace test

};  // namespace pace2024

#endif