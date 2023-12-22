#ifndef PACE2024_CROSSINGS_HPP
#define PACE2024_CROSSINGS_HPP

#include <vector>

#include "bipartite_graph.hpp"
#include "matrix.hpp"

namespace pace2024 {

template <typename T>
T compute_crossings(const general_bipartite_graph<T>& graph,
                    const folded_square_matrix<T>& cr_matrix,
                    const std::vector<T>& ordering) {
    const T n1 = graph.get_n1();
    assert(n1 == cr_matrix.get_m() && n1 == cr_matrix.get_n() && n1 == ordering.size());
    std::vector<T> positions(n1);

    // we compute the positions first such that
    // we are able to access the matrix cache friendly
    for (std::size_t i = 0; i < n1; ++i) {
        positions[ordering[i]] = i;
    }

    T nof_crossings = 0;
    for (std::size_t i = 0; i < n1; ++i) {
        for (std::size_t j = i + 1, j < n1; ++j) {
            if (positions[i] < positions[j])
                nof_crossings += cr_matrix(i, j);
            else
                nof_crossings += cr_matrix(j, i);
        }
    }
    return nof_crossings;
}

};  // namespace pace2024

#endif