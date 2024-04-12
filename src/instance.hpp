#ifndef PACE_INSTANCE_HPP
#define PACE_INSTANCE_HPP

#include <filesystem>
#include <limits>

#include "bipartite_graph.hpp"
#include "input.hpp"
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
    bipartite_graph<T> *graph_ptr{new bipartite_graph<T>()};

    /**
     * @brief the crossing number matrix
     */
    folded_matrix<R> *cr_matrix_ptr{nullptr};

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
    const fs::path filepath;

    /**
     * @brief constructs graph and crossing matrix based on input from `filepath`
     *
     * @param filepath filepath to instance
     */
    instance(const fs::path filepath) : filepath(filepath) {
        parse_input(filepath, *graph_ptr);
        cr_matrix_ptr = new folded_matrix<R>(graph_ptr->get_n_free());
        lower_bound = fill_crossing_matrix(*graph_ptr, *cr_matrix_ptr);
    }

    // delete copy and move constructor as well as copy and move assignment
    instance(const instance<T, R> &other) = delete;
    instance(instance<T, R> &&other) = delete;
    instance<T, R> &operator=(const instance<T, R> &other) = delete;
    instance<T, R> &operator=(instance<T, R> &&other) = delete;

    ~instance() {
        delete graph_ptr;
        graph_ptr = nullptr;
        delete cr_matrix_ptr;
        cr_matrix_ptr = nullptr;
    }

    const bipartite_graph<T> &graph() const { return *graph_ptr; }

    const folded_matrix<R> &cr_matrix() const { return *cr_matrix_ptr; }

    const uint32_t &get_lower_bound() const { return lower_bound; }
};

};  // namespace pace

#endif