#include "io/input.hpp"

//
// public methods
//

void pace::input::lift_ordering(const std::size_t i,
                                const std::vector<vertex_t> &subordering,
                                std::vector<vertex_t> &ordering) {
    assert(tried_split);
    assert(i < get_n_subgraphs());
    assert(ordering.size() == graph.get_n_free());

    const std::size_t begin = i > 0 ? borders[i - 1] : 0;
    const std::size_t end = borders[i];
    assert(end == subordering.size() + begin);

    for (std::size_t j = begin, k = 0; j < end; ++j, ++k) {
        ordering[j] = free_layer[begin + subordering[k]];
    }
}

bool pace::input::try_split() {
    if (tried_split) return get_n_subgraphs() >= 2;
    tried_split = true;
    if (graph.get_m() == 0 || graph.get_n_free() == 0) {
        first_graph_empty = true;
        return false;
    }
    sort_free_layer();

    // borders contains indices where instances begin
    const std::size_t n_free = graph.get_n_free();
    std::size_t i = n_free;
    --i;
    vertex_t leftmost_nbor = std::numeric_limits<vertex_t>::max();
    while (i > 0) {
        // since graph.get_m() > 0, and by the sorting, last vertex has degree >
        // 0
        leftmost_nbor =
            std::min(graph.get_leftmost_nbor(free_layer[i]), leftmost_nbor);

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

bool pace::input::compare(const vertex_t &u, const vertex_t &v) const {
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

//
// private methods
//

inline void pace::input::add_subgraph(  //
    const std::size_t begin, const std::size_t end) {
    assert(tried_split);
    assert(begin < end);

    // construct new bipartite graph and add to list
    bipartite_graph *subgraph_ptr = new bipartite_graph();
    subgraph_ptrs.emplace_back(subgraph_ptr);
    const std::size_t n_free = end - begin;
    subgraph_ptr->add_free_vertices(n_free);
    vertex_t v_free = 0;
    vertex_t v_fixed = 0;
    std::size_t m = 0;

    // setup map for fixed vertices
    std::map<vertex_t, vertex_t> map_fixed;
    for (std::size_t i = begin; i < end; ++i) {
        const vertex_t &v = free_layer[i];
        for (const vertex_t &u : graph.get_neighbors(v)) {
            map_fixed.insert({u, 0});
        }
    }
    for (auto &[v, u] : map_fixed) {
        u = v_fixed;
        ++v_fixed;
    }
    subgraph_ptr->add_fixed_vertices(v_fixed);

    // add edges to subgraph
    for (std::size_t i = begin; i < end; ++i) {
        const vertex_t &v = free_layer[i];
        for (const vertex_t &u : graph.get_neighbors(v)) {
            auto it = map_fixed.find(u);
            assert(it != map_fixed.end());
            const vertex_t u_fixed = it->second;
            subgraph_ptr->add_edge(u_fixed, v_free);
        }

        m += graph.get_degree(v);
        ++v_free;
    }
    subgraph_ptr->sort_adjacency_lists();

    assert(subgraph_ptr->is_sorted());
    assert(m == subgraph_ptr->get_m());
    assert(v_free == n_free);
    assert(map_fixed.size() == subgraph_ptr->get_n_fixed());
}

inline void pace::input::sort_free_layer() {
    const std::size_t n_free = graph.get_n_free();
    free_layer.resize(n_free);
    for (vertex_t v = 0; v < n_free; ++v) free_layer[v] = v;
    std::sort(free_layer.begin(), free_layer.end(),
              [=](const vertex_t &u, const vertex_t &v) -> bool {
                  return this->compare(u, v);
              });
}
