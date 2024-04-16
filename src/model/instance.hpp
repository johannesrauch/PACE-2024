#ifndef PACE_MODEL_INSTANCE_HPP
#define PACE_MODEL_INSTANCE_HPP

#include <limits>
#include <memory>

#include "matrix/matrix.hpp"
#include "model/bipartite_graph.hpp"

namespace pace {

namespace fs = std::filesystem;

/**
 * @brief class that comprises all relevant data for solving the instance
 */
class instance {
    /**
     * @brief the bipartite graph to the one-sided crossing minimization instance
     */
    const bipartite_graph &graph;

    /**
     * @brief the crossing number matrix (on demand)
     */
    std::unique_ptr<crossing_matrix> cr_matrix_ptr;

    /**
     * @brief lower bound to the optimal value of this instance
     */
    uint32_t lower_bound{0};

    /**
     * @brief upper bound to the optimal value of this instance
     */
    uint32_t upper_bound{std::numeric_limits<uint32_t>::max()};

   public:
    instance(const bipartite_graph &graph) : graph(graph) {}

    // delete copy constructor and assignment
    instance(const instance &other) = delete;
    instance &operator=(const instance &other) = delete;

    //
    // getter
    //

    const bipartite_graph &get_graph() const { return graph; }

    /**
     * @brief returns crossing number matrix (and constructs if non-existant)
     */
    const crossing_matrix &get_cr_matrix() {
        if (!cr_matrix_ptr) create_cr_matrix();
        return *cr_matrix_ptr;
    }

    /**
     * @brief returns const reference to lower_bound (and initializes it beforehand if necessary)
     *
     * @return const uint32_t&
     */
    const uint32_t &get_lower_bound() {
        if (!cr_matrix_ptr) create_cr_matrix();
        return lower_bound;
    }

    uint32_t &get_upper_bound() { return upper_bound; }

    //
    // private helper methods
    //

   private:
    void create_cr_matrix() {
        cr_matrix_ptr = std::make_unique<crossing_matrix>(graph.get_n_free());
        lower_bound = fill_crossing_matrix(graph, *cr_matrix_ptr);
    }
};

class instance_view {
   protected:
    instance &instance_;

    const bipartite_graph &graph;
    const std::size_t n_free;
    const std::size_t n_fixed;

    uint32_t &upper_bound;

   public:
    instance_view(instance &instance_)
        : instance_(instance_),
          graph(instance_.get_graph()),
          n_free(graph.get_n_free()),
          n_fixed(graph.get_n_fixed()),
          upper_bound(instance_.get_upper_bound()) {}

    const crossing_matrix &cr_matrix() { return instance_.get_cr_matrix(); }

    const uint32_t &lower_bound() { return instance_.get_lower_bound(); }
};

};  // namespace pace

#endif