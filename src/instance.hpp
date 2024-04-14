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

/**
 * @brief general input that is able to split into (smaller) instances
 *
 * @tparam T vertex type
 * @tparam R crossing number type
 */
template <typename T = uint16_t, typename R = uint32_t>
class input {
    bipartite_graph<T> *graph_ptr{new bipartite_graph<T>()};

    std::vector<instance<T, R> *> instances;

    std::size_t n_instances{1};

    bool trivial_instance{false};

   public:
    input(const fs::path filepath) {
        parse_input(filepath, *graph_ptr);

        // sort free_layer according to rightmost nbors
        const std::size_t n_free = graph_ptr->get_n_free();
        std::vector<T> free_layer(n_free);
        for (T v = 0; v < n_free; ++v) free_layer[v] = v;
        std::sort(free_layer.begin(), free_layer.end(), [=](const T &u, const T &v) -> bool {
            const bool deg_u_0 = graph_ptr->degree_of_free(u) == 0;
            const bool deg_v_0 = graph_ptr->degree_of_free(v) == 0;
            if (deg_u_0 && deg_v_0) return u < v;
            if (deg_u_0) return true;
            if (deg_v_0) return false;
            const T u_r = graph_ptr->get_rightmost_nbor_of_free(u);
            const T v_r = graph_ptr->get_rightmost_nbor_of_free(v);
            if (u_r < v_r) return true;
            if (u_r > v_r) return false;
            const T u_l = graph_ptr->get_leftmost_nbor_of_free(u);
            const T v_l = graph_ptr->get_leftmost_nbor_of_free(v);
            return u_l < v_l;
        });

        // leftmost_nbors[i] is the leftmost nbor of free_layer[i], ..., free_layer[n_free - 1]
        std::vector<T> leftmost_nbors(n_free);
        T leftmost = n_free;
        auto rit_free_layer = free_layer.rbegin();
        for (auto rit = leftmost_nbors.rbegin(); rit != leftmost_nbors.rend(); ++rit) {
            assert(rit_free_layer != free_layer.rend());
            if (graph_ptr->degree_of_free(*rit_free_layer) == 0) break;
            leftmost = std::min(leftmost, graph_ptr->get_leftmost_nbor_of_free(*rit_free_layer));
            *rit = leftmost;
            ++rit_free_layer;
        }

        // borders contains indices where subinstances in free_layer begin
        std::vector<T> borders;
        T i = 0;
        while (i < n_free && graph_ptr->degree_of_free(free_layer[i]) == 0) ++i;
        if (i > 0) {
            trivial_instance = true;
            borders.emplace_back(i);
        }
        while (i + 1u < n_free) {
            assert(graph_ptr->degree_of_free(free_layer[i]) > 0);
            const T rightmost_nbor = graph_ptr->get_rightmost_nbor_of_free(free_layer[i]);
            if (rightmost_nbor <= leftmost_nbors[i + 1u]) {
                borders.emplace_back(i + 1u);
            }
            ++i;
        }

        n_instances = borders.size() + 1u;
    }

    // delete copy and move constructor as well as copy and move assignment
    input(const input<T, R> &other) = delete;
    input(input<T, R> &&other) = delete;
    input<T, R> &operator=(const input<T, R> &other) = delete;
    input<T, R> &operator=(input<T, R> &&other) = delete;

    ~input() {
        for (const instance<T, R> *instance_ptr : instances) {
            delete instance_ptr;
        }
        instances.clear();
    }

    std::size_t get_n_instances() { return n_instances; }

    bool exists_trivial_instance() { return trivial_instance; }
};

};  // namespace pace

#endif