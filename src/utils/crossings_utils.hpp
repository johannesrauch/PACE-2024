#ifndef PACE_UTILS_CROSSINGS_UTILS_HPP
#define PACE_UTILS_CROSSINGS_UTILS_HPP

#include <algorithm>
#include <utility>
#include <vector>

#include "matrix/matrix.hpp"
#include "model/bipartite_graph.hpp"
#include "utils/vector_utils.hpp"

namespace pace {

/**
 * @brief returns the crossing numbers of the vertices u and v
 *
 * @tparam T vertex type of graph
 * @tparam R return type
 */
template <typename T, typename R = uint32_t>
std::pair<R, R> crossing_numbers_of(const general_bipartite_graph<T>& graph,
                                    T u, T v) {
    assert(graph.is_sorted());
    const auto& nbors_u = graph.get_neighbors(u);
    const auto& nbors_v = graph.get_neighbors(v);

    const std::size_t deg_u = nbors_u.size();
    const std::size_t deg_v = nbors_v.size();

    // cases with no crossings
    if (deg_u == 0 || deg_v == 0) return std::make_pair(0, 0);
    if (nbors_u[deg_u - 1] <= nbors_v[0])
        return std::pair<R, R>{
            0, deg_u * deg_v - (nbors_u[deg_u - 1] < nbors_v[0] ? 0 : 1)};
    if (nbors_v[deg_v - 1] <= nbors_u[0])
        return std::pair<R, R>{
            deg_u * deg_v - (nbors_v[deg_v - 1] < nbors_u[0] ? 0 : 1), 0};

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
template <typename T, typename R = crossing_number_t>
crossing_number_t number_of_crossings(  //
    const general_bipartite_graph<T>& graph, const std::vector<T>& ordering) {
    assert(graph.is_sorted());

    // compute the positions of each element (inverse of the permutation
    // ordering)
    const std::size_t n1 = graph.get_n_free();
    std::vector<T> positions(n1);
    inverse(ordering, positions);

    // sort, we need the positions of the ends of the edges in the free layer
    auto& edges = const_cast<general_bipartite_graph<T>&>(graph).get_edges();
    std::sort(edges.begin(), edges.end(),
              [&](const std::pair<T, T>& a, const std::pair<T, T>& b) {
                  if (a.first < b.first)
                      return true;
                  else if (a.first > b.first)
                      return false;
                  else
                      return positions[a.second] < positions[b.second];
              });

    // build the accumulator tree
    std::size_t q = ordering.size();
    std::size_t firstindex = 1;
    while (firstindex < q) firstindex *= 2;
    // number of tree nodes
    std::size_t treesize = 2 * firstindex - 1;
    // index of leftmost leaf
    firstindex -= 1;
    std::vector<R> tree(treesize);

    // count the crossings
    uint32_t n_crossings = 0;
    for (std::size_t k = 0; k < graph.get_m(); ++k) {
        // insert edge k
        std::size_t index = positions[edges[k].second] + firstindex;
        ++tree[index];
        while (index > 0) {
            if (index % 2) n_crossings += tree[index + 1];
            index = (index - 1) / 2;
            ++tree[index];
        }
    }

    return n_crossings;
}

/**
 * @brief returns the number of crossings in the given ordering
 *
 * @tparam T vertex type
 * @tparam R crossing number type
 * @param cr_matrix crossing number matrix
 * @param ordering ordering of free layer
 */
template <typename T, typename R>
uint32_t number_of_crossings(const folded_matrix<R>& cr_matrix,
                             const std::vector<T>& ordering) {
    const std::size_t n1 = cr_matrix.get_m();
    assert(n1 == ordering.size());
    std::vector<T> positions(n1);
    inverse(ordering, positions);  // to access `cr_matrix` cache-friendly

    uint32_t n_crossings = 0;
    for (std::size_t i = 0; i < n1; ++i) {
        for (std::size_t j = i + 1; j < n1; ++j) {
            if (positions[i] < positions[j])
                n_crossings += cr_matrix(i, j);
            else
                n_crossings += cr_matrix(j, i);
        }
    }
    return n_crossings;
}

/**
 * @brief returns sum min(c_uv, c_vu), which is a lower bound for the instance
 *
 * @tparam R crossing number type
 * @param cr_matrix crossing number matrix
 * @return uint32_t lower bound for the instance
 */
template <typename R>
uint32_t get_lower_bound(const folded_matrix<R>& cr_matrix) {
    const std::size_t n2 = cr_matrix.get_n2();
    uint32_t lb = 0;
    for (std::size_t i = 0; i < n2; i += 2) {
        lb += std::min(cr_matrix.get_element(i), cr_matrix.get_element(i + 1));
    }
    return lb;
}

};  // namespace pace

#endif