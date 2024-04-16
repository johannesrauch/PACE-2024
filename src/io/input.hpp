#ifndef PACE_IO_INPUT_HPP
#define PACE_IO_INPUT_HPP

#include <limits>
#include <unordered_map>
#include <vector>

#include "model/bipartite_graph.hpp"
#include "io/parse_input.hpp"

namespace pace {

/**
 * @brief general input that is able to split into (smaller) instances
 */
class input {
    bipartite_graph graph;
    std::vector<std::unique_ptr<bipartite_graph>> subgraph_ptrs;

    std::vector<vertex_t> free_layer;
    std::vector<vertex_t> borders;

    bool first_graph_empty{false};
    bool tried_split{false};

   public:
    const fs::path filepath;

    input(const fs::path filepath) : filepath(filepath) { parse_input(filepath, graph); }

    // delete copy constructor and assignment
    input(const input &other) = delete;
    input &operator=(const input &other) = delete;

    const bipartite_graph &get_graph() const { return graph; }

    std::size_t get_n_subgraphs() {
        if (!tried_split) try_split();
        return subgraph_ptrs.size();
    }

    const bipartite_graph &get_subgraph(const std::size_t i) {
        if (!tried_split) try_split();
        assert(i < get_n_subgraphs());
        return *subgraph_ptrs[i];
    }

    bool is_first_graph_empty() {
        if (!tried_split) try_split();
        return first_graph_empty;
    }

    bool try_split() {
        tried_split = true;
        if (graph.get_m() == 0) {
            first_graph_empty = true;
            return false;
        }
        sort_free_layer();

        // borders contains indices where instances begin
        borders.clear();
        const std::size_t n_free = graph.get_n_free();
        std::size_t i = n_free;
        if (i > 0) --i;
        vertex_t leftmost_nbor = std::numeric_limits<vertex_t>::max();
        while (i > 0) {
            // since graph.get_m() > 0, and by the sorting, last vertex has degree > 0
            leftmost_nbor = std::min(graph.get_leftmost_nbor(free_layer[i]), leftmost_nbor);

            --i;
            if (graph.get_degree(free_layer[i]) == 0) {
                borders.emplace_back(i + 1u);
                first_graph_empty = true;
                break;
            }

            const vertex_t rightmost_nbor = graph.get_rightmost_nbor(free_layer[i]);
            if (rightmost_nbor <= leftmost_nbor) {
                borders.emplace_back(i + 1u);
            }
        }
        std::reverse(borders.begin(), borders.end());
        borders.emplace_back(n_free);

        // return if we cannot split
        if (borders.size() == 1) return false;

        // otherwise construct the other graphs
        subgraph_ptrs.clear();
        i = 0;
        std::size_t begin = 0, end;
        while (i < borders.size()) {
            end = borders[i];
            add_subgraph(begin, end);
            begin = end;
            ++i;
        }
        return true;
    }

    bool compare(const vertex_t &u, const vertex_t &v) const {
        const bool deg_u_0 = graph.get_degree(u) == 0;
        const bool deg_v_0 = graph.get_degree(v) == 0;
        if (deg_u_0 && deg_v_0) return u < v;
        if (deg_u_0) return true;
        if (deg_v_0) return false;
        const vertex_t u_r = graph.get_rightmost_nbor(u);
        const vertex_t v_r = graph.get_rightmost_nbor(v);
        if (u_r < v_r) return true;
        if (u_r > v_r) return false;
        const vertex_t u_l = graph.get_leftmost_nbor(u);
        const vertex_t v_l = graph.get_leftmost_nbor(v);
        return u_l < v_l;
    }

   private:
    inline void add_subgraph(const std::size_t begin, const std::size_t end) {
        assert(begin < end);
        assert(tried_split);
        // construct new bipartite graph and add to list
        bipartite_graph *subgraph_ptr = new bipartite_graph();
        subgraph_ptrs.emplace_back(subgraph_ptr);

        // basic parameters
        const std::size_t n_free = end - begin;
        subgraph_ptr->add_free_vertices(n_free);
        vertex_t v_free = 0;

        std::unordered_map<vertex_t, vertex_t> map_fixed;
        for (std::size_t i = begin; i < end; ++i) {
            const vertex_t &v = free_layer[i];

            // add neighbors of v to subgraph
            const std::vector<vertex_t> &nbors_v = graph.get_neighbors(v);
            for (const vertex_t &u : nbors_v) {
                auto it = map_fixed.find(u);
                vertex_t u_fixed;
                if (it == map_fixed.end()) {
                    u_fixed = map_fixed.size();
                    map_fixed.insert({u, u_fixed});
                    subgraph_ptr->add_fixed_vertices(1);
                } else {
                    u_fixed = it->second;
                }
                subgraph_ptr->add_edge(u_fixed, v_free);
            }

            ++v_free;
        }
        assert(v_free == n_free);
        assert(map_fixed.size() == subgraph_ptr->get_n_fixed());
    }

    inline void sort_free_layer() {
        const std::size_t n_free = graph.get_n_free();
        free_layer.resize(n_free);
        for (vertex_t v = 0; v < n_free; ++v) free_layer[v] = v;
        std::sort(free_layer.begin(), free_layer.end(),
                  [=](const vertex_t &u, const vertex_t &v) -> bool { return this->compare(u, v); });
    }
};

};  // namespace pace

#endif
