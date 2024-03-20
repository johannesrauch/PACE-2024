#ifndef PACE_INSTANCE_HPP
#define PACE_INSTANCE_HPP

#include <filesystem>

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
    bipartite_graph<T> *graph_;
    folded_matrix<R> *cr_matrix_;

   public:
    instance(const fs::path filepath) {
        graph_ = new bipartite_graph<T>();
        parse_input(filepath, *graph_);
        cr_matrix_ = new folded_matrix<R>(graph_->get_n_free());
        fill_crossing_matrix(*graph_, *cr_matrix_);
    }

    // delete copy and move constructor as well as copy and move assignment
    instance(const instance<T, R> &other) = delete;
    instance(instance<T, R> &&other) = delete;
    instance<T, R> &operator=(const instance<T, R> &other) = delete;
    instance<T, R> &operator=(instance<T, R> &&other) = delete;

    ~instance() {
        delete graph_;
        delete cr_matrix_;
    }

    const bipartite_graph<T> &graph() const { return *graph_; }
    bipartite_graph<T> &graph() { return *graph_; }
    const folded_matrix<R> &cr_matrix() const { return *cr_matrix_; }
    folded_matrix<R> &cr_matrix() { return *cr_matrix_; }
};

};  // namespace pace

#endif