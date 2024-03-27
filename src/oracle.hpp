#ifndef PACE_ORACLE_HPP
#define PACE_ORACLE_HPP

#include <limits>

#include "crossings.hpp"
#include "digraph.hpp"
#include "index.hpp"
#include "instance.hpp"
#include "topological_sort.hpp"
#include "transitive_hull.hpp"
#include "vector_utils.hpp"

namespace pace {

enum pattern { indeterminate = 0, u_before_v, v_before_u };

/**
 * @brief oracle that tells you if we are able to fix u < v or v < u
 *
 * @tparam T vertex type
 * @tparam R crossing number type
 */
template <typename T, typename R>
class oracle {
    /// @brief instance graph
    const bipartite_graph<T> &graph;

    /// @brief crossing number matrix
    const folded_matrix<R> &cr_matrix;

    /// @brief lower and upper bound on the optimal value
    const uint32_t lower_bound, upper_bound;

    /// @brief n = graph.get_n_free()
    const std::size_t n, n_choose_2;

   public:
    /// @brief initializes the oracle
    oracle(const instance<T, R> &instance, const uint32_t upper_bound)
        : graph(instance.graph()),  //
          cr_matrix(instance.cr_matrix()),
          lower_bound(instance.get_lower_bound()),
          upper_bound(upper_bound),
          n(graph.get_n_free()),
          n_choose_2(n * (n - 1) / 2) {
        assert(n > 0);
        assert(lower_bound <= upper_bound);
    }

    // todo: vertices with equal neighborhoods

    /**
     * @brief builds the restriction (di)graph and the magic vector
     *
     * @param graph out parameter
     * @param magic out parameter
     */
    uint32_t build(digraph<T> &restriction_graph,  //
                   std::vector<int> &magic,        //
                   std::vector<std::pair<T, T>> &unsettled) {
        assert(restriction_graph.get_n() == n);

        // see if there are settlable pairs
        restriction_graph.clear_arcs();
        magic.clear();
        magic.resize(n_choose_2);
        uint32_t offset = 0;
        for (T u = 0; u < n; ++u) {
            for (T v = u + 1; v < n; ++v) {
                offset += fix_if_possible(restriction_graph, magic, u, v);
            }
        }

#ifndef NDEBUG
        // test
        std::vector<T> ordering;
        assert(topological_sort(restriction_graph, ordering));
#endif

        // transitive hull
        std::vector<std::pair<T, T>> new_arcs;
        pace::transitive_hull(restriction_graph, new_arcs);
        for (const auto &[u, v] : new_arcs) {
            restriction_graph.add_arc(u, v);
            if (u < v) {
                magic[flat_index(n, n_choose_2, u, v)] = -2;
            } else {
                magic[flat_index(n, n_choose_2, v, u)] = -1;
            }
            offset += cr_matrix(u, v);
        }

        // restriction graph is done, magic needs some last work
        restriction_graph.set_rollback_point();
        assert(topological_sort(restriction_graph, ordering));
        identify_unsettled(magic, unsettled);
        return offset;
    }

    /// @brief tells you if we are able to fix u < v or v < u
    pattern foresee(const T u, const T v) {
        assert(u < v);

        pattern p = based_on_degree(u, v);
        if (p != indeterminate) return p;

        const R &c_uv = cr_matrix(u, v);
        const R &c_vu = cr_matrix(v, u);
        if (c_uv == c_vu) return indeterminate;

        // c_uv != c_vu from her
        p = based_on_crossing_numbers(c_uv, c_vu);
        if (p != indeterminate) return p;

        p = based_on_bounds(c_uv, c_vu);
        if (p != indeterminate) return p;

        p = based_on_pattern(u, v, c_uv, c_vu);
        return p;
    }

    /// @brief if one of the degrees is zero, we may impose an arbitrary (but fixed) order
    inline pattern based_on_degree(const T &u, const T &v) {
        if (graph.degree_of_free(u) == 0) {
            return u_before_v;
        } else if (graph.degree_of_free(v) == 0) {
            return v_before_u;
        } else {
            return indeterminate;
        }
    }

    /**
     * @brief assumption: c_uv != c_vu.
     * - if c_uv == 0, we can fix u < v, and
     * - if c_vu == 0, we can fix v < u.
     */
    inline pattern based_on_crossing_numbers(const R &c_uv, const R &c_vu) {
        assert(c_uv != c_vu);
        if (c_uv == 0) return u_before_v;
        if (c_vu == 0) return v_before_u;
        return indeterminate;
    }

    /**
     * @brief if it were the other way around, we already would have more than upper_bound crossings
     */
    inline pattern based_on_bounds(const R &c_uv, const R &c_vu) {
        const uint32_t diff = upper_bound - lower_bound;
        if (c_uv > c_vu && c_uv - c_vu > diff) {
            return v_before_u;
        }
        if (c_vu > c_uv && c_vu - c_uv > diff) {
            return u_before_v;
        }
        return indeterminate;
    }

    /**
     * @brief based on an idea of Dujmovic et al, see https://doi.org/10.1016/j.jda.2006.12.008
     */
    inline pattern based_on_pattern(const T &u, const T &v, const R &c_uv, const R &c_vu) {
        const std::vector<T> &nbors_u = graph.get_neighbors_of_free(u);
        const std::vector<T> &nbors_v = graph.get_neighbors_of_free(v);
        // necessary condition for further statements
        if (nbors_u.size() != nbors_v.size()) return indeterminate;
        assert(nbors_u.size() > 0 && nbors_v.size() > 0);

        std::vector<T> nbors_uv;
        sorted_vector_union(nbors_u, nbors_v, nbors_uv);
        std::size_t i_u = 0, i_v = 0;
        R a = nbors_u.size(), b = nbors_v.size();
        for (std::size_t i = 0; i < nbors_uv.size(); ++i) {
            if (i_u < nbors_u.size() && nbors_uv[i] == nbors_u[i_u]) {
                assert(a > 0);
                --a;
            }
            if (i_v < nbors_v.size() && nbors_uv[i] == nbors_v[i_v]) {
                assert(b > 0);
                --b;
            }

            if (c_uv < c_vu && a > b) return indeterminate;
            if (c_uv > c_vu && a < b) return indeterminate;

            if (i_v < nbors_v.size() && nbors_uv[i] == nbors_v[i_v]) {
                ++a;
                ++i_v;
            }
            if (i_u < nbors_u.size() && nbors_uv[i] == nbors_u[i_u]) {
                ++b;
                ++i_u;
            }

            if (c_uv < c_vu && a > b) return indeterminate;
            if (c_uv > c_vu && a < b) return indeterminate;
        }
        assert(i_u == nbors_u.size());
        assert(i_v == nbors_v.size());
        assert(a == nbors_v.size());
        assert(b == nbors_u.size());

        if (c_uv < c_vu) {
            return u_before_v;
        } else {
            return v_before_u;
        }
    }

   private:
    inline uint32_t fix_if_possible(digraph<T> &restriction_graph, std::vector<int> &magic, const T &u, const T &v) {
        switch (foresee(u, v)) {
            case u_before_v: {
                restriction_graph.add_arc(u, v);
                magic[flat_index(n, n_choose_2, u, v)] = -2;
                return cr_matrix(u, v);
            }

            case v_before_u: {
                restriction_graph.add_arc(v, u);
                magic[flat_index(n, n_choose_2, u, v)] = -1;
                return cr_matrix(v, u);
            }

            default:
                return 0;
        }
    }

    inline void identify_unsettled(std::vector<int> &magic, std::vector<std::pair<T, T>> &unsettled) {
        std::size_t i = 0, j = 0;
        for (T u = 0; u < n; ++u) {
            for (T v = u + 1; v < n; ++v) {
                assert(i < n_choose_2);
                if (magic[i] >= 0) {
                    magic[i] = j;
                    unsettled.emplace_back(u, v);
                    ++j;
                }
                ++i;
            }
        }
        assert(i == n_choose_2);
    }
};

};  // namespace pace

#endif