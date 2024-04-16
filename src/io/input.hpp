#ifndef PACE_IO_INPUT_HPP
#define PACE_IO_INPUT_HPP

#include <limits>
#include <unordered_map>
#include <vector>

#include "bipartite_graph.hpp"
#include "io/parse_input.hpp"

namespace pace {

/**
 * @brief general input that is able to split into (smaller) instances
 *
 * @tparam T vertex type
 * @tparam R crossing number type
 */
template <typename T = uint16_t, typename R = uint32_t>
class input {
    bipartite_graph<T> graph;
    std::vector<std::unique_ptr<bipartite_graph<T>>> subgraph_ptrs;

    std::vector<T> free_layer;
    std::vector<T> borders;

    bool first_graph_empty{false};
    bool tried_split{false};

   public:
    const fs::path filepath;

    input(const fs::path filepath) : filepath(filepath) { parse_input(filepath, graph); }

    // delete copy constructor and assignment
    input(const input<T, R> &other) = delete;
    input<T, R> &operator=(const input<T, R> &other) = delete;

    bool compare(const T &u, const T &v) const {
        const bool deg_u_0 = graph.degree_of_free(u) == 0;
        const bool deg_v_0 = graph.degree_of_free(v) == 0;
        if (deg_u_0 && deg_v_0) return u < v;
        if (deg_u_0) return true;
        if (deg_v_0) return false;
        const T u_r = graph.get_rightmost_nbor(u);
        const T v_r = graph.get_rightmost_nbor(v);
        if (u_r < v_r) return true;
        if (u_r > v_r) return false;
        const T u_l = graph.get_leftmost_nbor(u);
        const T v_l = graph.get_leftmost_nbor(v);
        return u_l < v_l;
    }

    const bipartite_graph<T> &get_graph() const { return graph; }

    std::size_t get_n_subgraphs() {
        if (!tried_split) try_split();
        return subgraph_ptrs.size();
    }

    const bipartite_graph<T> &get_subgraph(const std::size_t i) {
        if (!tried_split) try_split();
        assert(i < get_n_subgraphs());
        return subgraph_ptrs;
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
        T leftmost_nbor = std::numeric_limits<T>::max();
        while (i > 0) {
            // since graph.get_m() > 0, and by the sorting, last vertex has degree > 0
            leftmost_nbor = std::min(graph.get_leftmost_nbor(free_layer[i]), leftmost_nbor);

            --i;
            if (graph.degree_of_free(free_layer[i]) == 0) {
                borders.emplace_back(i + 1u);
                first_graph_empty = true;
                break;
            }

            const T rightmost_nbor = graph.get_rightmost_nbor(free_layer[i]);
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

   private:
    inline void add_subgraph(const std::size_t begin, const std::size_t end) {
        assert(begin < end);
        assert(tried_split);
        // construct new bipartite graph and add to list
        bipartite_graph<T> *subgraph_ptr = new bipartite_graph<T>();
        subgraph_ptrs.emplace_back(subgraph_ptr);

        // basic parameters
        const std::size_t n_free = end - begin;
        subgraph_ptr->add_free_vertices(n_free);
        T v_free = 0;

        std::unordered_map<T, T> map_fixed;
        for (std::size_t i = begin; i < end; ++i) {
            const T &v = free_layer[i];

            // add neighbors of v to subgraph
            const std::vector<T> &nbors_v = graph.get_neighbors_of_free(v);
            for (const T &u : nbors_v) {
                auto it = map_fixed.find(u);
                T u_fixed;
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
        for (T v = 0; v < n_free; ++v) free_layer[v] = v;
        std::sort(free_layer.begin(), free_layer.end(),
                  [=](const T &u, const T &v) -> bool { return this->compare(u, v); });
    }
};

};  // namespace pace

#endif
