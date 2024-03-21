#ifndef PACE_BARYCENTER_HEURISTIC_HPP
#define PACE_BARYCENTER_HEURISTIC_HPP

#include <vector>

#include "crossings.hpp"
#include "instance.hpp"
#include "shift_heuristic.hpp"

namespace pace {

/**
 * @brief barycenter heuristic solver
 *
 * @tparam T vertey type
 * @tparam R crossing numbers type
 */
template <typename T, typename R>
class barycenter_heuristic {
    /// @brief input instance of one-sided crossing minimization
    const bipartite_graph<T> &graph;

    /// @brief crossing number matrix
    const folded_matrix<R> &cr_matrix;

    /// @brief number of vertices in the free layer
    const std::size_t n_free;

    /// @brief barycenters[u] = (sum over all neighbors v of u) / degree(u)
    std::vector<double> barycenters;

   public:
    /// @brief shift heuristic improver
    shift_heuristic<T, R> shift_h;

    /**
     * @brief sets internal references to graph and ordering
     *
     * @param graph input instance of one-sided crossing minimization
     * @param ordering in-out parameter, in which the computed ordering is stored
     */
    barycenter_heuristic(const instance<T, R> &instance)
        : graph(instance.graph()),
          cr_matrix(instance.cr_matrix()),
          n_free(graph.get_n_free()),
          barycenters(n_free),
          shift_h(instance) {
        fill_barycenters();
    }

    // delete copy and move constructor as well as copy and move assignment
    barycenter_heuristic(const barycenter_heuristic<T, R> &other) = delete;
    barycenter_heuristic(barycenter_heuristic<T, R> &&other) = delete;
    barycenter_heuristic<T, R> operator()(const barycenter_heuristic<T, R> &other) = delete;
    barycenter_heuristic<T, R> operator()(barycenter_heuristic<T, R> &&other) = delete;

    /// @brief runs the barycenter heuristic and stores the ordering in `ordering`
    template <bool SHIFT = true>
    uint32_t operator()(std::vector<T> &ordering) {
        ordering.resize(n_free);
        for (std::size_t i = 0; i < n_free; ++i) ordering[i] = i;
        std::sort(
            ordering.begin(),
            ordering.end(),
            [&](const T &a, const T &b) -> bool {
                if (barycenters[a] < barycenters[b])
                    return true;
                else if (barycenters[a] > barycenters[b])
                    return false;
                else
                    return a < b;
            });
        if constexpr (SHIFT) {
            return shift_h(ordering);
        } else {
            return number_of_crossings(graph, ordering);
        }
    }

   private:
    /// @brief initializes `ordering` and computes barycenters and stores them in `barycenters`
    inline void fill_barycenters() {
        for (T i = 0; i < n_free; ++i) barycenters[i] = barycenter(i);
    }

    /// @brief returns (sum of every neighbor j of i) / degree(i)
    inline double barycenter(const T &i) {
        const std::vector<T> &neighbors = graph.get_neighbors_of_free(i);
        if (neighbors.empty()) return 0.;
        double sum = 0.;
        for (const T &j : neighbors) sum += j;
        return sum / neighbors.size();
    }
};

};  // namespace pace

#endif