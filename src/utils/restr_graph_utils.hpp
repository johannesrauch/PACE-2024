#ifndef PACE_UTILS_RESTR_GRAPH_UTILS_HPP
#define PACE_UTILS_RESTR_GRAPH_UTILS_HPP

#include "utils/topological_sort.hpp"

namespace pace {

namespace internal {

template <typename T>
void build_restr_graph(const std::vector<double> &column_values,             //
                       const std::vector<std::pair<T, T>> &unsettled_pairs,  //
                       general_digraph<T> &restr_graph) {
    assert(column_values.size() == unsettled_pairs.size());
    restr_graph.rollback();
    std::size_t i = 0;
    for (const auto &[u, v] : unsettled_pairs) {
        const double x = column_values[i];
        if (x < 0.5) {
            restr_graph.add_arc(v, u);
        } else {
            restr_graph.add_arc(u, v);
        }
        ++i;
    }
}

};  // namespace internal

template <typename T>
void build_restr_graph_ordering(const std::vector<double> &column_values,             //
                             const std::vector<std::pair<T, T>> &unsettled_pairs,  //
                             general_digraph<T> &restr_graph,                      //
                             std::vector<T> &ordering) {
    internal::build_restr_graph(column_values, unsettled_pairs, restr_graph);
#ifndef NDEBUG
    bool acyclic =
#endif
        topological_sort(restr_graph, ordering);
    assert(acyclic);
}

};  // namespace pace

#endif
