#ifndef PACE2024_BIPARTITE_GRAPH_HPP
#define PACE2024_BIPARTITE_GRAPH_HPP

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "debug_printf.hpp"

namespace pace2024 {

/**
 * @brief a bipartite graph implementation
 * used for storing the pace instances
 *
 * @tparam T vertex type
 * @tparam std::enable_if_t<std::is_unsigned<T>::value>
 */
template <typename T, class = typename std::enable_if_t<std::is_unsigned<T>::value>>
class bipartite_graph {
   private:
    /**
     * @brief number of vertices in the fixed partite set
     */
    std::size_t n_fixed{0};

    /**
     * @brief number of vertices in the free partite set
     */
    std::size_t n_free{0};

    /**
     * @brief number of edges
     */
    std::size_t m{0};

    /**
     * @brief neighbors of vertices in the free layer; each adjacency list is sorted
     */
    std::vector<std::vector<T>> adjacency_lists;

    /**
     * @brief edges
     * first vertex is in fixed layer
     * second vertex is in free layer
     */
    std::vector<std::pair<T, T>> edges;

   public:
    using vertextype = T;

    bipartite_graph() {}

    // delete copy constructor, move constructor, copy assignment and move assignment
    bipartite_graph(const bipartite_graph &rhs) = delete;
    bipartite_graph(bipartite_graph &&rhs) = delete;
    bipartite_graph &operator=(const bipartite_graph &rhs) = delete;
    bipartite_graph &operator=(bipartite_graph &&rhs) = delete;

    /**
     * @brief adds the edge uv to the graph
     * 
     * @param u vertex of fixed layer
     * @param v vertex of free layer
     */
    void add_edge(const T u, const T v) {
        assert(u < n_fixed);
        assert(v < n_free);
        adjacency_lists[v].emplace_back(u);
        edges.emplace_back(u, v);
        ++m;
    }

    /**
     * @brief returns number of vertices in the fixed layer
     */
    std::size_t get_n_fixed() const { return n_fixed; }

    /**
     * @brief returns number of vertices in the free layer
     */
    std::size_t get_n_free() const { return n_free; }

    /**
     * @brief returns number of edges
     */
    std::size_t get_m() const { return m; }

    /**
     * @brief returns a constant reference to all adjacency lists
     */
    const std::vector<std::vector<T>> &get_adjacency_lists() const {
        return adjacency_lists;
    }

    /**
     * @brief get neighbors of vertex v (which is in the free layer)
     *
     * @param v
     * @return const std::vector<T>&
     */
    const std::vector<T> &get_neighbors(const T v) const {
        assert(v < n_free);
        return adjacency_lists[v];
    }

    /**
     * @brief returns a reference to all edges
     */
    std::vector<std::pair<T, T>> &get_edges() {
        return edges;
    }

    /**
     * @brief returns a constant reference to all edges
     */
    const std::vector<std::pair<T, T>> &get_edges() const {
        return edges;
    }

    /**
     * @brief sets the number of fixed vertices
     */
    void set_n_fixed(const std::size_t n_fixed_) {
        n_fixed = n_fixed_;
    }

    /**
     * @brief sets the number of free vertices; resizes adjacency list
     */
    void set_n_free(const std::size_t n_free_) {
        n_free = n_free_;
        adjacency_lists.resize(n_free);
    }

    /**
     * @brief sorts all adjacency lists in ascending order
     */
    void sort_adjacency_lists() {
        for (auto &adjacency_list : adjacency_lists) {
            std::sort(adjacency_list.begin(), adjacency_list.end());
        }
    }
};

using uint64_bipartite_graph = bipartite_graph<uint64_t>;
using uint32_bipartite_graph = bipartite_graph<uint32_t>;
using uint16_bipartite_graph = bipartite_graph<uint16_t>;

};  // namespace pace2024

#endif