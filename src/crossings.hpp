#ifndef PACE2024_CROSSINGS_HPP
#define PACE2024_CROSSINGS_HPP

#include <algorithm>
#include <utility>
#include <vector>

#include "bipartite_graph.hpp"
#include "matrix.hpp"
#include "vector_utils.hpp"

namespace pace2024 {

/**
 * @brief returns the crossing numbers of the vertices u and v
 *
 * @tparam T vertex type of graph
 * @tparam R return type
 */
template <typename T, typename R = uint32_t>
std::pair<R, R> crossing_numbers_of(const bipartite_graph<T>& graph, T u, T v) {
    const auto& nbors_u = graph.get_neighbors_of_free(u);
    const auto& nbors_v = graph.get_neighbors_of_free(v);

    const std::size_t deg_u = nbors_u.size();
    const std::size_t deg_v = nbors_v.size();

    // cases with no crossings
    if (deg_u == 0 || deg_v == 0) return std::make_pair(0, 0);
    if (nbors_u[deg_u - 1] <= nbors_v[0])
        return std::pair<R, R>{0, deg_u * deg_v - (nbors_u[deg_u - 1] < nbors_v[0] ? 0 : 1)};
    if (nbors_v[deg_v - 1] <= nbors_u[0])
        return std::pair<R, R>{deg_u * deg_v - (nbors_v[deg_v - 1] < nbors_u[0] ? 0 : 1), 0};

    // compute number of common neighbors of u and v
    R nof_common_nbors = sorted_vector_intersection<T, R>(nbors_u, nbors_v);

    // compute crossing number c_uv
    R c_uv = 0;
    for (std::size_t i = 0, j = 0; i < deg_u && j < deg_v;) {
        if (nbors_u[i] == nbors_v[j]) {
            ++i;
            ++j;
            c_uv += deg_u - i;
        } else if (nbors_u[i] < nbors_v[j]) {
            ++i;
        } else {
            ++j;
            c_uv += deg_u - i;
        }
    }

    assert(deg_u * deg_v >= nof_common_nbors + c_uv);
    return std::make_pair(c_uv, deg_u * deg_v - nof_common_nbors - c_uv);
}

/**
 * @brief computes the number of crossings of the given ordering.
 * from https://doi.org/10.1007/3-540-36151-0_13
 *
 * @tparam T vertex type
 * @tparam R accumulation tree type
 */
template <typename T, typename R = uint32_t>
uint32_t number_of_crossings(const bipartite_graph<T>& graph, const std::vector<T>& ordering) {
    // compute the positions of each element
    // (inverse of the permutation ordering)
    const std::size_t n1 = graph.get_n_free();
    std::vector<T> positions(n1);
    for (std::size_t i = 0; i < n1; ++i) {
        positions[ordering[i]] = i;
    }

    // sort, we need the positions of the ends of the edges in the free layer
    auto& edges = const_cast<bipartite_graph<T>&>(graph).get_edges();
    std::sort(edges.begin(), edges.end(), [&](const std::pair<T, T>& a, const std::pair<T, T>& b) {
        if (a.first < b.first)
            return true;
        else if (a.first > b.first)
            return false;
        else
            return positions[a.second] < positions[b.second];
    });

    /* build the accumulator tree */
    std::size_t q = ordering.size();
    std::size_t firstindex = 1;
    while (firstindex < q) firstindex *= 2;
    std::size_t treesize = 2 * firstindex - 1; /* number of tree nodes */
    firstindex -= 1;                           /* index of leftmost leaf */
    std::vector<R> tree(treesize);

    /* count the crossings */
    uint32_t crosscount = 0;                                 /* number of crossings */
    for (std::size_t k = 0; k < graph.get_m(); ++k) { /* insert edge k */
        std::size_t index = positions[edges[k].second] + firstindex;
        ++tree[index];
        while (index > 0) {
            if (index % 2) crosscount += tree[index + 1];
            index = (index - 1) / 2;
            ++tree[index];
        }
    }
    return crosscount;
}

/**
 * @brief returns the number of crossings in the given ordering
 *
 * @tparam T vertex type
 * @tparam R data type
 * @param cr_matrix crossing number matrix
 */
template <typename T, typename R>
uint32_t number_of_crossings(const folded_matrix<R>& cr_matrix, const std::vector<T>& ordering) {
    const std::size_t n1 = cr_matrix.get_m();
    assert(n1 == ordering.size());
    std::vector<T> positions(n1);
    inverse(ordering, positions);  // to access `cr_matrix` cache-friendly

    uint32_t nof_crossings = 0;
    for (std::size_t i = 0; i < n1; ++i) {
        for (std::size_t j = i + 1; j < n1; ++j) {
            if (positions[i] < positions[j])
                nof_crossings += cr_matrix(i, j);
            else
                nof_crossings += cr_matrix(j, i);
        }
    }
    return nof_crossings;
}

template <typename T, typename R>
uint32_t number_of_crossings(const sparse_matrix<R>& cr_matrix, const std::vector<T>& ordering) {
    const std::size_t n_free = cr_matrix.get_m();
    assert(n_free == ordering.size());
    std::vector<T> positions(n_free);
    inverse(ordering, positions);  // to access `cr_matrix` cache-friendly

    uint32_t nof_crossings = 0;
    for (std::size_t i = 0; i < cr_matrix.get_nof_nonzero_elements(); ++i) {
        const std::pair<T, T> &p = cr_matrix.get_indices(i);
        if (p.first < p.second && positions[p.first] < positions[p.second]) {
            nof_crossings += cr_matrix.get_datum(i);
        } else if (p.first > p.second && positions[p.first] < positions[p.second]) {
            nof_crossings += cr_matrix.get_datum(i);
        }
    }
    return nof_crossings;
}

};  // namespace pace2024

#endif