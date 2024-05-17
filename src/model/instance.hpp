#ifndef PACE_MODEL_INSTANCE_HPP
#define PACE_MODEL_INSTANCE_HPP

#include <chrono>
#include <limits>
#include <memory>

#include "log/debug_printf.hpp"
#include "matrix/matrix.hpp"
#include "model/bipartite_graph.hpp"
#include "model/digraph.hpp"
#include "utils/crossings_utils.hpp"
#include "utils/index_utils.hpp"
#include "utils/topological_sort.hpp"
#include "utils/transitive_hull.hpp"

namespace pace {

enum pattern { indeterminate = 0, u_before_v, v_before_u };

class highs_lp;

/**
 * @brief class that comprises all relevant data for solving the instance
 */
class instance {
   public:
    /**
     * @brief the bipartite graph to the one-sided crossing minimization
     * instance
     */
    const bipartite_graph &graph;
    const std::size_t n_free;
    const std::size_t n_free_2;
    const std::size_t n_fixed;

   private:
    /**
     * @brief the crossing number matrix (on demand)
     */
    std::unique_ptr<crossing_matrix> cr_matrix_unique_ptr;
    crossing_matrix *cr_matrix_ptr{nullptr};  // for subinstances

    /**
     * @brief arc u->v in restriction_graph implies u before v (on demand)
     */
    std::unique_ptr<digraph> restriction_graph_ptr;

    /**
     * @brief stores the problem kernel (on demand)
     * if magic[flat_index(u, v)] < 0, ~magic[flat_index(u, v)] is the value of
     * x_uv in an optimal solution. otherwise, magic[flat_index(u, v)] stores
     * the column index of x_uv in the lp.
     */
    std::vector<magic_t> magic;

    /**
     * @brief stores u < v s.t. magic[flat_index(u, v)] >= 0, that is, these
     * x_uv are in the lp (on demand)
     */
    std::vector<std::pair<vertex_t, vertex_t>> unsettled_pairs;

    /**
     * @brief objective value offset for the lp
     */
    crossing_number_t objective_offset{0};

    /**
     * @brief lower bound to the optimal value of this instance
     */
    crossing_number_t lower_bound{0};

    /**
     * @brief upper bound to the optimal value of this instance
     */
    crossing_number_t upper_bound{
        std::numeric_limits<crossing_number_t>::max()};

    /**
     * @brief ordering corresponding to upper_bound
     */
    std::vector<vertex_t> ordering;

    /**
     * @brief time when best solution was found
     */
    std::chrono::time_point<std::chrono::system_clock> t_sol{
        std::chrono::system_clock::now()};

   public:
    instance(const bipartite_graph &graph,
             const bool initialize_ordering = true)
        : graph(graph),
          n_free(graph.get_n_free()),
          n_free_2(n_free * (n_free - 1) / 2),
          n_fixed(graph.get_n_fixed()),
          ordering(n_free) {
        if (initialize_ordering) {
            identity(n_free, ordering);
            upper_bound = number_of_crossings(graph, ordering);
        }
        assert(n_free > 0);
    }

    // delete copy and move constructor and assignment
    instance(const instance &other) = delete;
    instance(instance &&other) = delete;
    instance &operator=(const instance &other) = delete;
    instance &operator=(instance &&other) = delete;

    //
    // getter
    //

    const bipartite_graph &get_graph() const { return graph; }

    const crossing_matrix &get_cr_matrix() {
        if (!cr_matrix_ptr) create_cr_matrix();
        return *cr_matrix_ptr;
    }

    digraph &get_restriction_graph() {
        if (!restriction_graph_ptr) create_kernel();
        return *restriction_graph_ptr;
    }

    const std::vector<magic_t> &get_magic() {
        if (!restriction_graph_ptr) create_kernel();
        return magic;
    }

    const std::vector<std::pair<vertex_t, vertex_t>> &get_unsettled_pairs() {
        if (!restriction_graph_ptr) create_kernel();
        return unsettled_pairs;
    }

    const crossing_number_t &get_objective_offset() {
        if (!restriction_graph_ptr) create_kernel();
        return objective_offset;
    }

    const crossing_number_t &get_lower_bound() {
        if (!cr_matrix_ptr) create_cr_matrix();
        return lower_bound;
    }

    const crossing_number_t &get_upper_bound() const { return upper_bound; }

    const std::vector<vertex_t> &get_ordering() const { return ordering; }

    const std::chrono::time_point<std::chrono::system_clock> &get_t_sol() {
        return t_sol;
    }

    //
    // setter
    //

    void update_ordering(const std::vector<vertex_t> &another_ordering,
                         const crossing_number_t n_crossings) {
        if (n_crossings < upper_bound) {
            upper_bound = n_crossings;
            ordering = another_ordering;
            t_sol = std::chrono::system_clock::now();
        }
    }

    void update_ordering(const std::vector<vertex_t> &another_ordering) {
        update_ordering(another_ordering,
                        number_of_crossings(graph, another_ordering));
    }

    void update_lower_bound(const crossing_number_t lb) {
        lower_bound = std::max(lower_bound, lb);
    }

    // void update_upper_bound(crossing_number_t ub) { upper_bound =
    // std::min(upper_bound, ub); }

    void update_kernel(const std::vector<double> &lp_sol);

    //
    // subinstance
    //

    instance *new_rins_instance(highs_lp &lp);
    instance *new_lsearch_instance(const uint16_t lsearch_width);

    //
    // private helper methods
    //

   private:
    instance *new_subinstance();

    void create_cr_matrix() {
        cr_matrix_unique_ptr =
            std::make_unique<crossing_matrix>(graph.get_n_free());
        cr_matrix_ptr = cr_matrix_unique_ptr.get();
        lower_bound = fill_crossing_matrix(graph, *cr_matrix_ptr);
    }

    void create_kernel(highs_lp *lp = nullptr,
                       const uint16_t lsearch_width = 0);

    void create_unsettled_pairs();

    inline void fix_u_before_v(const vertex_t &u, const vertex_t &v) {
        restriction_graph_ptr->add_arc(u, v);
        magic[flat_index(n_free, n_free_2, u, v)] = -2;
        objective_offset += (*cr_matrix_ptr)(u, v);
    }

    inline void fix_v_before_u(const vertex_t &u, const vertex_t &v) {
        restriction_graph_ptr->add_arc(v, u);
        magic[flat_index(n_free, n_free_2, u, v)] = -1;
        objective_offset += (*cr_matrix_ptr)(v, u);
    }

    void compute_transitive_hull();

    //
    // problem kernel methods
    //

    /**
     * @brief tells if we are able to fix u < v or v < u
     */
    pattern foresee(const vertex_t &u, const vertex_t &v, highs_lp *lp,
                    const std::vector<vertex_t> &positions,
                    const uint16_t &lsearch_width);

    /**
     * @brief if we construct the kernel for a rins instance, we fix were the
     * relaxation agrees with our current best solution
     */
    pattern based_on_relaxation(const vertex_t &u, const vertex_t &v,
                                highs_lp &lp,
                                const std::vector<vertex_t> &positions);

    pattern based_on_ordering(const vertex_t &u, const vertex_t &v,
                              const std::vector<vertex_t> &positions,
                              const uint16_t &lsearch_width);

   public:
    /**
     * @brief if one of the degrees is zero, we may impose an arbitrary (but
     * fixed) order
     */
    pattern based_on_degree(const vertex_t &u, const vertex_t &v);

    /**
     * @brief assumption: c_uv != c_vu.
     * - if c_uv == 0, we can fix u < v, and
     * - if c_vu == 0, we can fix v < u.
     */
    pattern based_on_crossing_numbers(const crossing_number_t &c_uv,
                                      const crossing_number_t &c_vu);

    /**
     * @brief if it were the other way around, we already would have more than
     * upper_bound crossings
     */
    pattern based_on_bounds(const crossing_number_t &c_uv,
                            const crossing_number_t &c_vu);

    /**
     * @brief based on an idea of Dujmovic et al, see
     * https://doi.org/10.1016/j.jda.2006.12.008
     */
    pattern based_on_pattern(const vertex_t &u,              //
                             const vertex_t &v,              //
                             const crossing_number_t &c_uv,  //
                             const crossing_number_t &c_vu);
};

struct instance_view {
    instance &instance_;

    const bipartite_graph &graph;
    const std::size_t &n_free;
    const std::size_t &n_free_2; /**< = n_free choose 2 */
    const std::size_t &n_fixed;

    const crossing_number_t &upper_bound;

    instance_view(instance &instance_)
        : instance_(instance_),
          graph(instance_.get_graph()),
          n_free(instance_.n_free),
          n_free_2(instance_.n_free_2),
          n_fixed(instance_.n_fixed),
          upper_bound(instance_.get_upper_bound()) {}

    const crossing_matrix &cr_matrix() { return instance_.get_cr_matrix(); }

    const crossing_number_t &lower_bound() {
        return instance_.get_lower_bound();
    }

    digraph &restriction_graph() { return instance_.get_restriction_graph(); }

    const std::vector<std::pair<vertex_t, vertex_t>> &unsettled_pairs() {
        return instance_.get_unsettled_pairs();
    }

    const std::vector<magic_t> &magic() { return instance_.get_magic(); }

    const crossing_number_t &objective_offset() {
        return instance_.get_objective_offset();
    }

    const std::vector<vertex_t> &get_ordering() const {
        return instance_.get_ordering();
    }

    void update_ordering(const std::vector<vertex_t> &another_ordering) {
        instance_.update_ordering(another_ordering);
    }

    void update_ordering(const std::vector<vertex_t> &another_ordering,
                         const crossing_number_t n_crossings) {
        instance_.update_ordering(another_ordering, n_crossings);
    }

    void update_lower_bound(const crossing_number_t lb) {
        instance_.update_lower_bound(lb);
    }

    // void update_upper_bound(crossing_number_t ub) {
    // instance_.update_upper_bound(ub); }
};

};  // namespace pace

#endif
