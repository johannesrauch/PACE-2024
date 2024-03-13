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
template <typename T, class = typename std::enable_if_t<std::is_integral<T>::value>>
class bipartite_graph {
   private:
    /**
     * @brief neighbors of vertices in the free layer
     */
    std::vector<std::vector<T>> adjacency_lists_of_fixed;

    /**
     * @brief neighbors of vertices in the free layer
     */
    std::vector<std::vector<T>> adjacency_lists_of_free;

    /**
     * @brief edges
     * first vertex is in fixed layer
     * second vertex is in free layer
     */
    std::vector<std::pair<T, T>> edges;

   public:
    using vertextype = T;

    bipartite_graph() {}

    bipartite_graph(const bipartite_graph &rhs)
        : adjacency_lists_of_fixed(rhs.adjacency_lists_of_fixed),
          adjacency_lists_of_free(rhs.adjacency_lists_of_free),
          edges(rhs.edges) {
        // assert that copy constructor is only called on empty graphs to avoid copying
        assert(rhs.get_n() == 0);
    }

    // delete move constructor, copy assignment and move assignment
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
        adjacency_lists_of_fixed[u].emplace_back(v);
        adjacency_lists_of_free[v].emplace_back(u);
        edges.emplace_back(u, v);
    }

    /**
     * @brief adds a vertex to the fixed layer
     */
    std::size_t add_fixed_vertex() {
        const std::size_t ret = get_n_fixed();
        adjacency_lists_of_fixed.emplace_back();
        return ret;
    }

    /**
     * @brief adds `n` vertices to the fixed layer
     */
    void add_fixed_vertices(const std::size_t n) {
        adjacency_lists_of_fixed.resize(adjacency_lists_of_fixed.size() + n);
    }

    /**
     * @brief adds a vertex to the free layer
     */
    std::size_t add_free_vertex() {
        const std::size_t ret = get_n_free();
        adjacency_lists_of_free.emplace_back();
        return ret;
    }

    /**
     * @brief adds `n` vertices to the free layer
     *
     * @param n
     */
    void add_free_vertices(const std::size_t n) {
        adjacency_lists_of_free.resize(adjacency_lists_of_free.size() + n);
    }

    //
    // clear modifiers
    //

    /**
     * @brief clears and resets all attributes of this objects
     */
    void clear() {
        adjacency_lists_of_fixed.clear();
        adjacency_lists_of_free.clear();
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
    std::size_t get_n_fixed() const { return adjacency_lists_of_fixed.size(); }

    /**
     * @brief returns number of vertices in the free layer
     */
    std::size_t get_n_free() const { return adjacency_lists_of_free.size(); }

    /**
     * @brief returns number of edges
     */
    std::size_t get_m() const { return edges.size(); }

    /**
     * @brief get neighbors of vertex v, which is in the fixed layer
     */
    const std::vector<T> &get_neighbors_of_fixed(const T v) const {
        assert(v < get_n_fixed());
        return adjacency_lists_of_fixed[v];
    }

    /**
     * @brief get neighbors of vertex v, which is in the free layer
     */
    const std::vector<T> &get_neighbors_of_free(const T v) const {
        assert(v < get_n_free());
        return adjacency_lists_of_free[v];
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

    //
    // sorting methods
    //

    /**
     * @brief sorts all adjacency lists in ascending order
     */
    void sort_adjacency_lists() {
        for (auto &adjacency_list : adjacency_lists_of_fixed) {
            std::sort(adjacency_list.begin(), adjacency_list.end());
        }
        for (auto &adjacency_list : adjacency_lists_of_free) {
            std::sort(adjacency_list.begin(), adjacency_list.end());
        }
    }
};

using uint64_bipartite_graph = bipartite_graph<uint64_t>;
using uint32_bipartite_graph = bipartite_graph<uint32_t>;
using uint16_bipartite_graph = bipartite_graph<uint16_t>;

};  // namespace pace2024

#endif