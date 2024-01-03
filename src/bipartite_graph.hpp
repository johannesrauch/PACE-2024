#ifndef PACE2024_BIPARTITE_GRAPH_HPP
#define PACE2024_BIPARTITE_GRAPH_HPP

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <vector>

#include "instance.hpp"

namespace pace2024 {

template <typename T, class = typename std::enable_if_t<std::is_unsigned<T>::value>>
class general_bipartite_graph {
   private:
    std::size_t n0,                               // number of vertices in A, the fixed partite set
        n1,                                       // number of vertices in B, the free partite set
        m;                                        // number of edges
    std::vector<std::vector<T>> adjacency_lists;  // each adjacency list is sorted

   public:
    using datatype = T;

    general_bipartite_graph(const general_bipartite_graph &rhs) = delete;
    general_bipartite_graph &operator=(const general_bipartite_graph &rhs) = delete;
    general_bipartite_graph &operator=(general_bipartite_graph &&rhs) = delete;

    template <typename T1>
    general_bipartite_graph(const general_instance<T1> &instance)
        : n0(instance.get_n0()),
          n1(instance.get_n1()),
          m(instance.get_m()),
          adjacency_lists(n1) {
        auto edges = instance.get_edges();

        for (auto [x, y] : edges) {
            --x;  // normalize to 0, ..., n0-1
            assert(x < n0);
            y -= n0 + 1;  // normalize to 0, ..., n1-1
            assert(y < n1);
            adjacency_lists[y].emplace_back(x);
        }

        sort_adjacency_lists();
    }

    std::size_t get_n0() const { return n0; }

    std::size_t get_n1() const { return n1; }

    std::size_t get_m() const { return m; }

    const std::vector<std::vector<T>> &get_adjacency_lists() const {
        return adjacency_lists;
    }

   private:
    void sort_adjacency_lists() {
        for (auto &adjacency_list : adjacency_lists) {
            std::sort(adjacency_list.begin(), adjacency_list.end());
        }
    }
};

using uint64_bipartite_graph = general_bipartite_graph<std::uint64_t>;
using uint32_bipartite_graph = general_bipartite_graph<std::uint32_t>;

};  // namespace pace2024

#endif