#ifndef PACE_INSTANCE_HPP
#define PACE_INSTANCE_HPP

#include <filesystem>
#include <limits>
#include <memory>
#include <unordered_map>

#include "bipartite_graph.hpp"
#include "io/parse_input.hpp"
#include "matrix.hpp"

namespace pace {

namespace fs = std::filesystem;

/**
 * @brief class that comprises all relevant data for solving the instance
 *
 * @tparam T vertex type
 * @tparam R crossing numbers type
 */
template <typename T = uint16_t, typename R = uint32_t>
class instance {
    /**
     * @brief the bipartite graph to the one-sided crossing minimization instance
     */
    const bipartite_graph<T> &graph;

    /**
     * @brief the crossing number matrix
     */
    std::unique_ptr<folded_matrix<R>> cr_matrix_ptr;

    /**
     * @brief lower bound to the optimal value of this instance
     */
    uint32_t lower_bound{0};

    /**
     * @brief upper bound to the optimal value of this instance
     */
    uint32_t upper_bound{std::numeric_limits<uint32_t>::max()};

    /**
     * @brief an ordering of the vertices of graph
     */
    std::vector<T> ordering;

   public:
    instance(const bipartite_graph<T> &graph) : graph(graph) {}

    // delete copy constructor and assignment
    instance(const instance<T, R> &other) = delete;
    instance<T, R> &operator=(const instance<T, R> &other) = delete;

    const bipartite_graph<T> &get_graph() const { return graph; }

    const folded_matrix<R> &get_cr_matrix() const {
        if (!cr_matrix_ptr) {
            cr_matrix_ptr = std::make_unique<folded_matrix<R>>(graph);
        }
        return *cr_matrix_ptr;
    }

    const uint32_t &get_lower_bound() const { return lower_bound; }
};

};  // namespace pace

#endif