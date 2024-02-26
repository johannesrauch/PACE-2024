#ifndef PACE2024_GRAPH_HPP
#define PACE2024_GRAPH_HPP

#include <cstdint>
#include <queue>
#include <vector>

#include "debug.hpp"

namespace pace2024 {

template <typename T>
class general_graph {
   private:
    std::vector<std::vector<T>> adjacency_lists;
    std::vector<std::vector<T>> in_adjacency_lists;

   public:
    using datatype = T;

    general_graph(const general_graph &rhs) = delete;
    general_graph &operator=(const general_graph &rhs) = delete;
    general_graph &operator=(general_graph &&rhs) = delete;

    general_graph(const std::size_t n) : adjacency_lists(n), in_adjacency_lists(n) {}

    /**
     * @brief add arc (u,v)
     *
     * @param u
     * @param v
     */
    void add_edge(const T u, const T v) {
        assert(u < get_n() && v < get_n());
        adjacency_lists[u].emplace_back(v);
        in_adjacency_lists[v].emplace_back(u);
    }

    void pop_neighbor(const T u) {
        assert(u < get_n());
        if (adjacency_lists[u].size() > 0) adjacency_lists[u].pop_back();
    }

    void pop_in_neighbor(const T v) {
        assert(v < get_n());
        if (in_adjacency_lists[v].size() > 0) in_adjacency_lists[v].pop_back();
    }

    std::size_t get_n() const {
        return adjacency_lists.size();
    }

    const std::vector<T> &get_adjacency_list(const T v) const {
        return adjacency_lists[v];
    }

    const std::vector<T> &get_in_adjacency_list(const T v) const {
        return in_adjacency_lists[v];
    }

    /**
     * @brief checks if u and v are connected by a directed path
     *
     * @param u vertex, 0 <= u < get_n()
     * @param v vertex, 0 <= v < get_n()
     * @return true if u and v are connected by a directed path
     * @return false if not
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
 * 
 * @tparam T vertex type
 * @param graph 
 * @return true if the graph is a tournament
 * @return false if not
 */
template <typename T>
bool is_tournament(const general_graph<T> &graph) {
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