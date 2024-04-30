#ifndef PACE_UTILS_RESTR_GRAPH_UTILS_HPP
#define PACE_UTILS_RESTR_GRAPH_UTILS_HPP

#include <cmath>

#include "exact/highs_lp.hpp"
#include "utils/topological_sort.hpp"

namespace pace {

namespace internal {

/**
 * @brief see https://math.stackexchange.com/a/4166138.
 */
inline double smooth_transition(const double x, const double k = 2.) {
    assert(x > 1e-6);
    assert(x < 1. - 1e-6);
    const double y = 0.5 * (std::tanh(k * (2 * x - 1) / (2 * std::sqrt(x * (1 - x)))) + 1);
    assert(y >= 0);
    assert(y <= 1);
    return y;
}

};  // namespace internal

template <typename T>
void build_restr_graph(const std::vector<double> &column_values,             //
                       const std::vector<std::pair<T, T>> &unsettled_pairs,  //
                       general_digraph<T> &restr_graph) {
    assert(column_values.size() == unsettled_pairs.size());
    restr_graph.rollback();
    std::size_t i = 0;
    for (const auto &[u, v] : unsettled_pairs) {
        const double x = column_values[i];
        const bool add_uv = x > 1e-6 && (x >= 1. - 1e-6 || coinflip(internal::smooth_transition(x)));
        if (add_uv)
            restr_graph.add_arc(u, v);
        else
            restr_graph.add_arc(v, u);
        ++i;
    }
}

template <typename T>
bool build_restr_graph_ordering(const std::vector<double> &column_values,             //
                                const std::vector<std::pair<T, T>> &unsettled_pairs,  //
                                general_digraph<T> &restr_graph,                      //
                                std::vector<T> &ordering) {
    build_restr_graph(column_values, unsettled_pairs, restr_graph);
    return topological_sort(restr_graph, ordering);
}

};  // namespace pace

#endif
