#ifndef PACE_PATTERNS_HPP
#define PACE_PATTERNS_HPP

#include "bipartite_graph.hpp"
#include "matrix.hpp"
#include "vector_utils.hpp"

namespace pace {

enum pattern { indeterminate = 0, u_before_v, v_before_u };

template <typename T, typename R>
class pattern_analyzer {
    const bipartite_graph<T> &graph;

    const folded_matrix<R> &cr_matrix;

   public:
    pattern_analyzer(const bipartite_graph<T> &graph,  //
                     const folded_matrix<R> &cr_matrix)
        : graph(graph),  //
          cr_matrix(cr_matrix) {}

    pattern analyze(const T u, const T v) {
        const R &c_uv = cr_matrix(u, v);
        const R &c_vu = cr_matrix(v, u);
        if (c_uv == 0 && c_vu != 0) return u_before_v;
        if (c_uv != 0 && c_vu == 0) return v_before_u;
        if (c_uv == c_vu) return indeterminate;

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

        if (c_uv < c_vu)
            return u_before_v;
        else
            return v_before_u;
    }
};

};  // namespace pace

#endif