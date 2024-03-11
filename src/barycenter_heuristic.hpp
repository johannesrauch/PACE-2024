#ifndef PACE2024_BARYCENTER_HEURISTIC_HPP
#define PACE2024_BARYCENTER_HEURISTIC_HPP

#include <vector>

#include "bipartite_graph.hpp"

namespace pace2024 {

/**
 * @brief heuristic solver class that implements the barycenter heuristic with `run()`
 *
 * @tparam T vertey type
 */
template <typename T>
class barycenter_heuristic {
    /// @brief input instance of one-sided crossing minimization
    const bipartite_graph<T> &graph;

    /// @brief number of vertices in the free layer
    const std::size_t n1;

    /// @brief in-out vector for storing the ordering
    std::vector<T> &ordering;

    /// @brief barycenters[u] = (sum over all neighbors v of u) / degree(u)
    std::vector<double> barycenters;

   public:
    /**
     * @brief sets internal references to graph and ordering
     *
     * @param graph input instance of one-sided crossing minimization
     * @param ordering in-out parameter, in which the computed ordering is stored
     */
    barycenter_heuristic(const bipartite_graph<T> &graph,
                         std::vector<T> &ordering)
        : graph(graph),
          n1(graph.get_n_free()),
          ordering(ordering) {
    }

    // delete copy and move constructor as well as copy and move assignment
    barycenter_heuristic(const barycenter_heuristic<T> &other) = delete;
    barycenter_heuristic(barycenter_heuristic<T> &&other) = delete;
    barycenter_heuristic<T> operator()(const barycenter_heuristic<T> &other) = delete;
    barycenter_heuristic<T> operator()(barycenter_heuristic<T> &&other) = delete;

    /// @brief runs the barycenter heuristic and stores the ordering in `ordering`
    void run() {
        initialize();
        std::sort(ordering.begin(),
                  ordering.end(),
                  [&](const T &a, const T &b) -> bool {
                      if (barycenters[a] < barycenters[b])
                          return true;
                      else if (barycenters[a] > barycenters[b])
                          return false;
                      else
                          return a < b;
                  });
    }

   private:
    /// @brief initializes `ordering` and computes barycenters and stores them in `barycenters`
    inline void initialize() {
        ordering.resize(n1);
        for (T i = 0; i < n1; ++i) ordering[i] = i;
        barycenters.resize(n1);
        for (T i = 0; i < n1; ++i) barycenters[i] = barycenter(i);
    }

    /// @brief returns (sum of every neighbor j of i) / degree(i)
    inline double barycenter(const T &i) {
        const std::vector<T> &neighbors = graph.get_neighbors(i);
        if (neighbors.empty()) return 0.;
        double sum = 0.;
        for (const T &j : neighbors) sum += j;
        return sum / neighbors.size();
    }
};

};  // namespace pace2024

#endif