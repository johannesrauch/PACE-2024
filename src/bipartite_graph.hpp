#ifndef PACE2024_BIPARTITE_GRAPH_HPP
#define PACE2024_BIPARTITE_GRAPH_HPP

#include <algorithm>
#include <cstdint>
#include <vector>

#include "instance.hpp"

namespace pace2024 {

class bipartite_graph {
   private:
    uint64_t n0,  // number of vertices in A, the fixed partite set
        n1,       // number of vertices in B
        m;        // number of edges
    std::vector<std::vector<uint64_t>> adjacency_lists;

   public:
    bipartite_graph(const bipartite_graph &rhs) = delete;
    bipartite_graph &operator=(const bipartite_graph &rhs) = delete;
    bipartite_graph &operator=(bipartite_graph &&rhs) = delete;

    bipartite_graph(const instance &instance)
        : n0(instance.get_n0()),
          n1(instance.get_n1()),
          m(instance.get_m()),
          adjacency_lists(n1) {
        auto edges = instance.get_edges();

        for (auto [x, y] : edges) {
            --x;          // normalize to 0, ..., n0-1
            y -= n0 + 1;  // normalize to 0, ..., n1-1
            adjacency_lists[y].emplace_back(x);
        }
    }

    uint64_t get_n0() const { return n0; }

    uint64_t get_n1() const { return n1; }

    uint64_t get_m() const { return m; }

    const std::vector<std::vector<uint64_t>> &get_adjacency_lists() const {
        return adjacency_lists;
    }

    void sort_adjacency_lists() {
        for (auto &adjacency_list : adjacency_lists) {
            std::sort(adjacency_list.begin(), adjacency_list.end());
        }
    }
};

};  // namespace pace2024

#endif