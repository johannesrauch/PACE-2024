#ifndef PACE_MODEL_INSTANCE_HPP
#define PACE_MODEL_INSTANCE_HPP

#include <limits>
#include <memory>

#include "log/debug_printf.hpp"
#include "matrix/matrix.hpp"
#include "model/bipartite_graph.hpp"
#include "model/digraph.hpp"
#include "utils/index_utils.hpp"
#include "utils/topological_sort.hpp"
#include "utils/transitive_hull.hpp"

namespace pace {

enum pattern { indeterminate = 0, u_before_v, v_before_u };

/**
 * @brief class that comprises all relevant data for solving the instance
 */
class instance {
   public:
    /**
     * @brief the bipartite graph to the one-sided crossing minimization instance
     */
    const bipartite_graph &graph;
    const std::size_t n_free;
    const std::size_t n_free_2;
    const std::size_t n_fixed;

   private:
    /**
     * @brief the crossing number matrix (on demand)
     */
    std::unique_ptr<crossing_matrix> cr_matrix_ptr;

    /**
     * @brief arc u->v in restriction_graph implies u before v (on demand)
     */
    std::unique_ptr<digraph> restriction_graph_ptr;

    /**
     * @brief stores the problem kernel (on demand)
     * if magic[flat_index(u, v)] < 0, ~magic[flat_index(u, v)] is the value of x_uv in an optimal
     * solution. otherwise, magic[flat_index(u, v)] stores the column index of x_uv in the lp.
     */
    std::vector<magic_t> magic;

    /**
     * @brief stores u < v s.t. magic[flat_index(u, v)] >= 0, that is, these x_uv are in the lp (on
     * demand)
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
    crossing_number_t upper_bound{std::numeric_limits<uint32_t>::max()};

   public:
    instance(const bipartite_graph &graph)
        : graph(graph),
          n_free(graph.get_n_free()),
          n_free_2(n_free * (n_free - 1) / 2),
          n_fixed(graph.get_n_fixed()) {
        assert(n_free > 0);
    }

    // delete copy constructor and assignment
    instance(const instance &other) = delete;
    instance &operator=(const instance &other) = delete;

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

    const crossing_number_t &get_upper_bound() { return upper_bound; }

    //
    // setter
    //

    void update_upper_bound(crossing_number_t ub) { upper_bound = std::min(upper_bound, ub); }

    //
    // private helper methods
    //

   private:
    void create_cr_matrix() {
        cr_matrix_ptr = std::make_unique<crossing_matrix>(graph.get_n_free());
        lower_bound = fill_crossing_matrix(graph, *cr_matrix_ptr);
    }

    void create_kernel() {
        if (restriction_graph_ptr) return;

        PACE_DEBUG_PRINTF("start create_kernel\n");
        get_cr_matrix();  // initializes crossing matrix
        const std::size_t n_free = graph.get_n_free();
        restriction_graph_ptr = std::make_unique<digraph>(graph.get_n_free());
        magic.resize(n_free_2);

        // let's foresee if we can fix u < v or v < u
        for (vertex_t u = 0u; u + 1u < n_free; ++u) {
            for (vertex_t v = u + 1u; v < n_free; ++v) {
                pattern p = foresee(u, v);
                if (p == u_before_v) {
                    restriction_graph_ptr->add_arc(u, v);
                    magic[flat_index(n_free, n_free_2, u, v)] = -2;
                    objective_offset += (*cr_matrix_ptr)(u, v);
                } else if (p == v_before_u) {
                    restriction_graph_ptr->add_arc(v, u);
                    magic[flat_index(n_free, n_free_2, u, v)] = -1;
                    objective_offset += (*cr_matrix_ptr)(v, u);
                }
            }
        }

#ifndef NDEBUG
        std::vector<vertex_t> ordering;
        assert(topological_sort(*restriction_graph_ptr, ordering));
#endif

        // transitive hull
        std::vector<std::pair<vertex_t, vertex_t>> new_arcs;
        pace::transitive_hull(*restriction_graph_ptr, new_arcs);
        for (const auto &[u, v] : new_arcs) {
            restriction_graph_ptr->add_arc(u, v);
            if (u < v) {
                magic[flat_index(n_free, n_free_2, u, v)] = -2;
            } else {
                magic[flat_index(n_free, n_free_2, v, u)] = -1;
            }
            objective_offset += (*cr_matrix_ptr)(u, v);
        }

        // restriction graph is done, magic needs some last work
        restriction_graph_ptr->set_rollback_point();
        assert(topological_sort(*restriction_graph_ptr, ordering));
        std::size_t i = 0, j = 0;
        for (vertex_t u = 0u; u + 1u < n_free; ++u) {
            for (vertex_t v = u + 1u; v < n_free; ++v) {
                assert(i < n_free_2);
                if (magic[i] >= 0) {
                    magic[i] = j;
                    unsettled_pairs.emplace_back(u, v);
                    ++j;
                }
                ++i;
            }
        }
        assert(i == n_free_2);
        PACE_DEBUG_PRINTF("end   create_kernel\n");
    }

    //
    // problem kernel methods
    //

    /**
     * @brief tells if we are able to fix u < v or v < u
     */
    pattern foresee(const vertex_t &u, const vertex_t &v) {
        assert(u < v);

        pattern p = based_on_degree(u, v);
        if (p != indeterminate) return p;

        const crossing_number_t &c_uv = (*cr_matrix_ptr)(u, v);
        const crossing_number_t &c_vu = (*cr_matrix_ptr)(v, u);
        if (c_uv == c_vu) return pattern::indeterminate;

        // c_uv != c_vu from her
        p = based_on_crossing_numbers(c_uv, c_vu);
        if (p != indeterminate) return p;

        p = based_on_bounds(c_uv, c_vu);
        if (p != indeterminate) return p;

        p = based_on_pattern(u, v, c_uv, c_vu);
        return p;
    }

   public:
    /**
     * @brief if one of the degrees is zero, we may impose an arbitrary (but fixed) order
     */
    inline pattern based_on_degree(const vertex_t &u, const vertex_t &v) {
        if (graph.get_degree(u) == 0)
            return pattern::u_before_v;
        else if (graph.get_degree(v) == 0)
            return pattern::v_before_u;
        else
            return pattern::indeterminate;
    }

    /**
     * @brief assumption: c_uv != c_vu.
     * - if c_uv == 0, we can fix u < v, and
     * - if c_vu == 0, we can fix v < u.
     */
    inline pattern based_on_crossing_numbers(const crossing_number_t &c_uv,
                                             const crossing_number_t &c_vu) {
        assert(c_uv != c_vu);
        if (c_uv == 0) return pattern::u_before_v;
        if (c_vu == 0) return pattern::v_before_u;
        return pattern::indeterminate;
    }

    /**
     * @brief if it were the other way around, we already would have more than upper_bound crossings
     */
    inline pattern based_on_bounds(const crossing_number_t &c_uv, const crossing_number_t &c_vu) {
        const crossing_number_t diff = upper_bound - lower_bound;
        if (c_uv > c_vu && c_uv - c_vu > diff) return pattern::v_before_u;
        if (c_vu > c_uv && c_vu - c_uv > diff) return pattern::u_before_v;
        return pattern::indeterminate;
    }

    /**
     * @brief based on an idea of Dujmovic et al, see https://doi.org/10.1016/j.jda.2006.12.008
     */
    inline pattern based_on_pattern(const vertex_t &u,              //
                                    const vertex_t &v,              //
                                    const crossing_number_t &c_uv,  //
                                    const crossing_number_t &c_vu) {
        const std::vector<vertex_t> &nbors_u = graph.get_neighbors(u);
        const std::vector<vertex_t> &nbors_v = graph.get_neighbors(v);
        // necessary condition for further statements
        if (nbors_u.size() != nbors_v.size()) return pattern::indeterminate;
        assert(nbors_u.size() > 0 && nbors_v.size() > 0);

        std::vector<vertex_t> nbors_uv;
        sorted_vector_union(nbors_u, nbors_v, nbors_uv);
        std::size_t i_u = 0, i_v = 0;
        crossing_number_t a = nbors_u.size(), b = nbors_v.size();

        for (std::size_t i = 0; i < nbors_uv.size(); ++i) {
            if (i_u < nbors_u.size() && nbors_uv[i] == nbors_u[i_u]) {
                assert(a > 0);
                --a;
            }
            if (i_v < nbors_v.size() && nbors_uv[i] == nbors_v[i_v]) {
                assert(b > 0);
                --b;
            }

            if (c_uv < c_vu && a > b) return pattern::indeterminate;
            if (c_uv > c_vu && a < b) return pattern::indeterminate;

            if (i_v < nbors_v.size() && nbors_uv[i] == nbors_v[i_v]) {
                ++a;
                ++i_v;
            }
            if (i_u < nbors_u.size() && nbors_uv[i] == nbors_u[i_u]) {
                ++b;
                ++i_u;
            }

            if (c_uv < c_vu && a > b) return pattern::indeterminate;
            if (c_uv > c_vu && a < b) return pattern::indeterminate;
        }

        assert(i_u == nbors_u.size());
        assert(i_v == nbors_v.size());
        assert(a == nbors_v.size());
        assert(b == nbors_u.size());

        if (c_uv < c_vu)
            return pattern::u_before_v;
        else
            return pattern::v_before_u;
    }
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

    const crossing_number_t &lower_bound() { return instance_.get_lower_bound(); }

    digraph &restriction_graph() { return instance_.get_restriction_graph(); }

    const std::vector<std::pair<vertex_t, vertex_t>> &unsettled_pairs() {
        return instance_.get_unsettled_pairs();
    }

    const std::vector<magic_t> &magic() { return instance_.get_magic(); }

    const crossing_number_t &objective_offset() { return instance_.get_objective_offset(); }

    void update_upper_bound(crossing_number_t ub) { instance_.update_upper_bound(ub); }
};

};  // namespace pace

#endif