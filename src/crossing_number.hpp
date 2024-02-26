#ifndef PACE2024_CROSSING_NUMBER_HPP
#define PACE2024_CROSSING_NUMBER_HPP

#include <algorithm>
#include <utility>
#include <vector>

#include "bipartite_graph.hpp"
#include "matrix.hpp"
#include "vector_intersection.hpp"

namespace pace2024 {

/**
 * @brief computes and returns the crossing numbers of the vertices u and v
 *
 * @tparam T vertex type of graph
 * @tparam R return type
 * @param graph
 * @param u vertex
 * @param v vertex
 * @return std::pair<R, R> crossing numbers of u and v
 */
template <typename T, typename R>
std::pair<R, R> crossing_numbers_of(const general_bipartite_graph<T>& graph, T u, T v) {
    auto& adjacency_lists = graph.get_adjacency_lists();
    auto& nbors_u = adjacency_lists[u];
    auto& nbors_v = adjacency_lists[v];

    const std::size_t deg_u = nbors_u.size();
    const std::size_t deg_v = nbors_v.size();

    // cases with no crossings
    if (deg_u == 0 || deg_v == 0) return std::make_pair(0, 0);
    if (nbors_u[deg_u - 1] <= nbors_v[0])
        return std::pair<R, R>{0, deg_u * deg_v - (nbors_u[deg_u - 1] < nbors_v[0] ? 0 : 1)};
    if (nbors_v[deg_v - 1] <= nbors_u[0])
        return std::pair<R, R>{deg_u * deg_v - (nbors_v[deg_v - 1] < nbors_u[0] ? 0 : 1), 0};

    R nof_common_nbors = vector_intersection<T, R>(nbors_u, nbors_v);

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

    return std::make_pair(c_uv, deg_u * deg_v - nof_common_nbors - c_uv);
}

/**
 * @brief computes the number of crossings of the given ordering
 *
 * based on https://doi.org/10.1007/3-540-36151-0_13
 *
 * @tparam T vertex type
 * @tparam R return type
 * @param ordering
 * @return R
 */
template <typename T, typename R>
R crossing_number_of(const general_bipartite_graph<T>& graph, const std::vector<T>& ordering) {
    // compute the positions of each element
    // (inverse of the permutation ordering)
    std::size_t n1 = graph.get_n1();
    std::vector<T> positions(n1);
    for (std::size_t i = 0; i < n1; ++i) {
        positions[ordering[i]] = i;
    }

    // sort for southside
    std::size_t m = graph.get_m();
    auto& edges = graph.get_edges();
    std::vector<T> southsequence(m);
    for (std::size_t i = 0; i < m; ++i) {
        southsequence[i] = i;
    }
    std::sort(southsequence.begin(),
              southsequence.end(),
              [&](const T& a, const T& b) {
                  if (edges[a].first < edges[b].first)
                      return true;
                  else if (edges[a].first < edges[b].first)
                      return false;
                  else
                      return positions[edges[a].second] < positions[edges[b].second];
              });

    /* build the accumulator tree */
    std::size_t q = ordering.size();
    std::size_t firstindex = 1;
    while (firstindex < q) firstindex *= 2;
    std::size_t treesize = 2 * firstindex - 1; /* number of tree nodes */
    firstindex -= 1;                           /* index of leftmost leaf */
    std::vector<R> tree(treesize);

    /* count the crossings */
    R crosscount = 0;                                 /* number of crossings */
    for (std::size_t k = 0; k < graph.get_m(); ++k) { /* insert edge k */
        std::size_t index = southsequence[edges[southsequence[k]].second] + firstindex;
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
 * @tparam T
 * @tparam R
 * @param cr_matrix crossing number matrix
 * @param ordering ordering of the free layer
 * @return R
 */
template <typename T, typename R>
R crossing_number_of(const folded_square_matrix<R>& cr_matrix, const std::vector<T>& ordering) {
    const std::size_t n1 = cr_matrix.get_m();
    assert(n1 == ordering.size());
    std::vector<T> positions(n1);

    // compute the positions first
    // then are able to access the matrix cache friendly
    for (std::size_t i = 0; i < n1; ++i) {
        positions[ordering[i]] = i;
    }

    R nof_crossings = 0;
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

};  // namespace pace2024

#endif