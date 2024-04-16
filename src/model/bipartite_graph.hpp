#ifndef PACE_MODEL_BIPARTITE_GRAPH_HPP
#define PACE_MODEL_BIPARTITE_GRAPH_HPP

#include <algorithm>
#include <cassert>
#include <type_traits>
#include <vector>

#include "types/types.hpp"

namespace pace {

/**
 * @brief a bipartite graph implementation
 * used for storing the pace instances
 *
 * @tparam T vertex type
 * @tparam std::enable_if_t<std::is_unsigned<T>::value>
 */
template <typename T, class = typename std::enable_if_t<std::is_unsigned<T>::value>>
class general_bipartite_graph {
    /**
     * @brief number of vertices in the fixed layer
     */
    std::size_t n_fixed{0};

    /**
     * @brief neighbors of vertices in the free layer
     */
    std::vector<std::vector<T>> adjacency_lists;

    /**
     * @brief edges, first vertex is in fixed layer, second vertex is in free layer
     */
    std::vector<std::pair<T, T>> edges;

    /**
     * @brief true iff all adjacency lists are sorted
     */
    bool is_sorted{false};

   public:
    using vertex_t = T;

    general_bipartite_graph() {}

    // delete copy and move constructor, copy assignment and move assignment
    general_bipartite_graph(const general_bipartite_graph &rhs) = delete;
    general_bipartite_graph(general_bipartite_graph &&rhs) = delete;
    general_bipartite_graph &operator=(const general_bipartite_graph &rhs) = delete;
    general_bipartite_graph &operator=(general_bipartite_graph &&rhs) = delete;

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
        is_sorted = false;
        adjacency_lists[v].emplace_back(u);
        edges.emplace_back(u, v);
    }

    /**
     * @brief adds `n` vertices to the fixed layer
     */
    void add_fixed_vertices(const std::size_t n) { n_fixed += n; }

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
    const std::vector<T> &get_neighbors(const T v) const {
        assert(v < get_n_free());
        return adjacency_lists[v];
    }

    /**
     * @brief returns the first element in neighbor list of v
     */
    T get_leftmost_nbor(const T v) const {
        assert(is_sorted);
        const std::vector<T> &nbors = get_neighbors(v);
        assert(nbors.size() > 0);
        return *get_neighbors(v).begin();
    }

    /**
     * @brief returns the last element in neighbor list of v
     */
    T get_rightmost_nbor(const T v) const {
        assert(is_sorted);
        const std::vector<T> &nbors = get_neighbors(v);
        assert(nbors.size() > 0);
        return *get_neighbors(v).rbegin();
    }

    /**
     * @brief returns a reference to all edges
     */
    std::vector<std::pair<T, T>> &get_edges() { return edges; }

    /**
     * @brief returns the degree of vertex `v` in the free layer
     */
    std::size_t get_degree(const std::size_t v) const { return get_neighbors(v).size(); }

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
        is_sorted = true;
    }
};

using bipartite_graph = general_bipartite_graph<vertex_t>;

};  // namespace pace

#endif