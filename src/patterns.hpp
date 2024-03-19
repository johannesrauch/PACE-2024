#ifndef PACE2024_PATTERNS_HPP
#define PACE2024_PATTERNS_HPP

#include "bipartite_graph.hpp"
#include "vector_utils.hpp"

namespace pace2024 {

enum pattern { indeterminate = 0, u_before_v, v_before_u };

template <typename T>
void analyze_pattern(const bipartite_graph<T> &graph, const T u, const T v) {
    const std::vector<T> &nbors_u = graph.get_neighbors_of_free(u);
    const std::vector<T> &nbors_v = graph.get_neighbors_of_free(v);
    if (nbors_u.size() == 0 || nbors_v.size() == 0) return;
    std::vector<T> nbors_uv;
    sorted_vector_union(nbors_u, nbors_v, nbors_uv);

    std::vector<T> a, b;
    a.reserve(2 * nbors_uv.size());
    b.reserve(2 * nbors_uv.size());
    std::size_t i = 0, i_u = 0, i_v = 0;
}

};  // namespace pace2024

#endif