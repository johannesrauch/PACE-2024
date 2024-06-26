#ifndef PACE_HEURISTIC_BARYCENTER_HEURISTIC_HPP
#define PACE_HEURISTIC_BARYCENTER_HEURISTIC_HPP

#include "heuristic/shift_heuristic.hpp"

namespace pace {

/**
 * @brief barycenter heuristic solver
 */
class barycenter_heuristic : public instance_view {
    /**
     * @brief barycenters[u] = (sum over all neighbors v of u) / degree(u)
     */
    std::vector<double> barycenters;
    shift_heuristic shift_h;

   public:
    /**
     * @brief sets internal references to graph and ordering
     *
     * @param graph input instance of one-sided crossing minimization
     * @param ordering in-out parameter, in which the computed ordering is stored
     */
    barycenter_heuristic(instance &instance_)
        : instance_view(instance_), barycenters(n_free), shift_h(instance_) {
        fill_barycenters();
    }

    // delete copy and move constructor as well as copy and move assignment
    barycenter_heuristic(const barycenter_heuristic &other) = delete;
    barycenter_heuristic(barycenter_heuristic &&other) = delete;
    barycenter_heuristic operator()(const barycenter_heuristic &other) = delete;
    barycenter_heuristic operator()(barycenter_heuristic &&other) = delete;

    /**
     * @brief runs the barycenter heuristic and stores the ordering in `ordering`
     */
    crossing_number_t operator()(std::vector<vertex_t> &ordering);

    inline bool compare(const vertex_t &a, const vertex_t &b) {
        if (barycenters[a] < barycenters[b])
            return true;
        else if (barycenters[a] > barycenters[b])
            return false;
        else
            return a < b;
    }

   private:
    /**
     * @brief initializes `ordering` and computes barycenters and stores them in `barycenters`
     */
    inline void fill_barycenters() {
        for (vertex_t i = 0; i < n_free; ++i) barycenters[i] = barycenter(i);
    }

    /**
     * @brief returns (sum of every neighbor j of i) / degree(i)
     */
    inline double barycenter(const vertex_t &i) {
        const std::vector<vertex_t> &neighbors = graph.get_neighbors(i);
        if (neighbors.empty()) return 0.;
        double sum = 0.;
        for (const vertex_t &j : neighbors) sum += j;
        return sum / neighbors.size();
    }
};

};  // namespace pace

#endif