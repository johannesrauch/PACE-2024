#ifndef PACE_DIGRAPH_HPP
#define PACE_DIGRAPH_HPP

#include <cstdint>
#include <vector>

#include "debug_printf.hpp"

namespace pace {

/**
 * @brief class for a generic directed graph with
 * definable vertex type
 *
 * @tparam T vertex type
 */
template <typename T = uint16_t>
class digraph {
   private:
    /// @brief stores all neighbors of each vertex
    std::vector<std::vector<T>> adjacency_lists;

   public:
    using vertex_type = T;

    // delete copy constructor, move constructor, copy assignment and move assignment
    digraph(const digraph &rhs) = delete;
    digraph(digraph &&rhs) = delete;
    digraph &operator=(const digraph &rhs) = delete;
    digraph &operator=(digraph &&rhs) = delete;

    /**
     * @brief constructs an empty digraph with n vertices
     *
     * @param n number of vertices
     */
    digraph(const std::size_t n) : adjacency_lists(n) {}

    /**
     * @brief add arc (u, v)
     *
     * @param u vertex, 0 <= u < get_n()
     * @param v vertex, 0 <= v < get_n()
     */
    void add_arc(const T u, const T v) {
        assert(u != v);
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
    const std::vector<T> &get_neighbors(const T v) const {
        return adjacency_lists[v];
    }

    /**
     * @brief deletes all arcs
     */
    void clear_arcs() {
        for (std::size_t i = 0; i < get_n(); ++i) {
            adjacency_lists[i].clear();
        }
    }
};

};  // namespace pace

#endif