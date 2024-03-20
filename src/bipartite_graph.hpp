#ifndef PACE_BIPARTITE_GRAPH_HPP
#define PACE_BIPARTITE_GRAPH_HPP

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "debug_printf.hpp"

namespace pace {

/**
 * @brief a bipartite graph implementation
 * used for storing the pace instances
 *
 * @tparam T vertex type
 * @tparam std::enable_if_t<std::is_unsigned<T>::value>
 */
template <typename T = uint16_t, class = typename std::enable_if_t<std::is_unsigned<T>::value>>
class bipartite_graph {
    /**
     * @brief number of vertices in the fixed layer
     */
    std::size_t n_fixed{0};

    /**
     * @brief neighbors of vertices in the free layer
     */
    std::vector<std::vector<T>> adjacency_lists;

    /**
     * @brief edges
     * first vertex is in fixed layer
     * second vertex is in free layer
     */
    std::vector<std::pair<T, T>> edges;

   public:
    using vertex_t = T;

    bipartite_graph() {}

    // delete copy and move constructor, copy assignment and move assignment
    bipartite_graph(const bipartite_graph &rhs) = delete;
    bipartite_graph(bipartite_graph &&rhs) = delete;
    bipartite_graph &operator=(const bipartite_graph &rhs) = delete;
    bipartite_graph &operator=(bipartite_graph &&rhs) = delete;

    //
    // add modifiers
    //

    /**
     * @brief adds the edge uv to the graph
     *
     * @param u vertex of fixed layer
     * @param v vertex of free layer
     */
    void add_edge(const T u, const T v) {
        assert(u < get_n_fixed());
        assert(v < get_n_free());
        adjacency_lists[v].emplace_back(u);
        edges.emplace_back(u, v);
    }

    /**
     * @brief adds `n` vertices to the fixed layer
     */
    void add_fixed_vertices(const std::size_t n) {
        n_fixed += n;
    }

    /**
     * @brief adds `n` vertices to the free layer
     */
    void add_free_vertices(const std::size_t n) {
        adjacency_lists.resize(adjacency_lists.size() + n);
    }

    //
    // clear modifiers
    //

    /**
     * @brief clears and resets all attributes of this objects
     */
    void clear() {
        adjacency_lists.clear();
        edges.clear();
    }

    //
    // getter
    //

    /**
     * @brief returns the number of all vertices
     */
    std::size_t get_n() const { return get_n_fixed() + get_n_free(); }

    /**
     * @brief returns number of vertices in the fixed layer
     */
    std::size_t get_n_fixed() const { return n_fixed; }

    /**
     * @brief returns number of vertices in the free layer
     */
    std::size_t get_n_free() const { return adjacency_lists.size(); }

    /**
     * @brief returns number of edges
     */
    std::size_t get_m() const { return edges.size(); }

    /**
     * @brief get neighbors of vertex v, which is in the free layer
     */
    const std::vector<T> &get_neighbors_of_free(const T v) const {
        assert(v < get_n_free());
        return adjacency_lists[v];
    }

    /**
     * @brief returns a reference to all edges
     */
    std::vector<std::pair<T, T>> &get_edges() { return edges; }

    /**
     * @brief returns a constant reference to all edges
     */
    const std::vector<std::pair<T, T>> &get_edges() const { return edges; }

    /**
     * @brief returns the degree of vertex `v` in the free layer
     */
    std::size_t degree_of_free(const std::size_t v) const {
        return get_neighbors_of_free(v).size();
    }

    //
    // sorting methods
    //

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

};  // namespace pace

#endif