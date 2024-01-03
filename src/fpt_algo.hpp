#ifndef PACE2024_FPT_ALGO_HPP
#define PACE2024_FPT_ALGO_HPP

#include <cassert>
#include <set>
#include <vector>

#include "bipartite_graph.hpp"
#include "graph.hpp"
#include "instance.hpp"
#include "matrix.hpp"
#include "median_heuristic.hpp"

namespace pace2024 {

namespace internal {

template <typename T>
struct pair_with_crossings {
    const T i, j, c_ij, c_ji, c;

    pair_with_crossings(T i, T j, T c_ij, T c_ji)
        : i(i),
          j(j),
          c_ij(c_ij),
          c_ji(c_ji),
          c(c_ij + c_ji) {
        assert(i < j);
    }
};

template <typename T>
struct pair_with_crossings_comparator {
    bool operator()(const pair_with_crossings<T> &a, const pair_with_crossings<T> &b) const {
        return a.c < b.c;
    }
};

};  // namespace internal

template <typename T>
class fpt_algo {
   private:
    const std::size_t n1;
    const general_bipartite_graph<T> bipartite_graph;
    const folded_square_matrix<T> cr_matrix;

    std::vector<T> ordering;
    T upper_bound;

    std::set<internal::pair_with_crossings<T>, internal::pair_with_crossings_comparator<T>> unsettled_pairs;
    pace2024::general_graph<T> acyclic_graph;

   public:
    fpt_algo(const general_instance<T> &instance)
        : n1(instance.get_n1()),
          bipartite_graph(instance),
          cr_matrix(bipartite_graph),
          ordering(n1),
          upper_bound(prob_median_heuristic(bipartite_graph, cr_matrix, ordering)),
          unsettled_pairs(),
          acyclic_graph(n1) {
        for (T i = 0; i < n1 - 1; ++i) {
            for (T j = i + 1; j < n1; ++j) {
                T c_ij = cr_matrix(i, j);
                T c_ji = cr_matrix(j, i);

                if (c_ij == 0) {
                    acyclic_graph.add_edge(i, j);  // fix i < j
                } else if (c_ji == 0) {
                    acyclic_graph.add_edge(j, i);  // fix j < i
                } else {
                    unsettled_pairs.emplace(i, j, c_ij, c_ji);
                }
            }
        }

        // todo: implement other data reduction rules
    }

    fpt_algo(const fpt_algo &rhs) = delete;
    fpt_algo &operator=(const fpt_algo &rhs) = delete;
    fpt_algo &operator=(fpt_algo &&rhs) = delete;
};

};  // namespace pace2024

#endif