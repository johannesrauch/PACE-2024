#include "model/instance.hpp"

#include "exact/highs_lp.hpp"
#include "utils/index_utils.hpp"

namespace pace {

//
// public methods
//

void instance::update_kernel(const std::vector<double> &lp_sol) {
    if (!restriction_graph_ptr) create_kernel();

    std::vector<vertex_t> positions;
    inverse(ordering, positions);

    const double tol = 1e-6;
    for (const auto &[u, v] : unsettled_pairs) {
        const std::size_t i = flat_index(n_free, n_free_2, u, v);
        if (positions[u] < positions[v] && lp_sol[i] >= 1. - tol) {
            fix_u_before_v(u, v);
        } else if (positions[v] < positions[u] && lp_sol[i] <= tol) {
            fix_v_before_u(u, v);
        }
    }

#ifndef NDEBUG
    std::vector<vertex_t> ordering;
    assert(topological_sort(*restriction_graph_ptr, ordering));
#endif

    compute_transitive_hull();

    // restriction graph is finished
    restriction_graph_ptr->set_rollback_point();
    assert(topological_sort(*restriction_graph_ptr, ordering));

    // magic and unsettled_pairs needs some last work
    create_unsettled_pairs();
}

//
// private methods
//

void pace::instance::create_kernel(highs_lp *lp) {
    if (restriction_graph_ptr) return;

    PACE_DEBUG_PRINTF("start create_kernel\n");
    get_cr_matrix();  // initializes crossing matrix
    const std::size_t n_free = graph.get_n_free();
    restriction_graph_ptr = std::make_unique<digraph>(graph.get_n_free());
    magic.resize(n_free_2);
    std::vector<vertex_t> positions;
    if (lp != nullptr) inverse(ordering, positions);

    // let's foresee if we can fix u < v or v < u
    for (vertex_t u = 0u; u + 1u < n_free; ++u) {
        for (vertex_t v = u + 1u; v < n_free; ++v) {
            switch (foresee(u, v, lp, positions)) {
                case u_before_v:
                    fix_u_before_v(u, v);
                    break;
                case v_before_u:
                    fix_v_before_u(u, v);
                    break;
                case indeterminate:
                    break;
            }
        }
    }

#ifndef NDEBUG
    std::vector<vertex_t> ordering;
    assert(topological_sort(*restriction_graph_ptr, ordering));
#endif

    compute_transitive_hull();

    // restriction graph is finished
    restriction_graph_ptr->set_rollback_point();
    assert(topological_sort(*restriction_graph_ptr, ordering));

    // magic and unsettled_pairs needs some last work
    create_unsettled_pairs();
    PACE_DEBUG_PRINTF("end   create_kernel\n");
}

void instance::create_unsettled_pairs() {
    assert(magic.size() == n_free_2);
    unsettled_pairs.clear();
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
}

void instance::compute_transitive_hull() {
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
}

pattern instance::foresee(const vertex_t &u, const vertex_t &v, highs_lp *lp,
                          const std::vector<vertex_t> &positions) {
    assert(u < v);

    pattern p = indeterminate;
    if (lp != nullptr) {
        p = based_on_relaxation(u, v, *lp, positions);
        if (p != indeterminate) return p;
    }

    p = based_on_degree(u, v);
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

pattern instance::based_on_relaxation(const vertex_t &u, const vertex_t &v,
                                      highs_lp &lp,
                                      const std::vector<vertex_t> &positions) {
    const std::size_t uv = flat_index(n_free, n_free_2, u, v);
    if (!lp.is_variable_integral(uv)) return indeterminate;
    const double x = lp.get_variable_value(uv);
    if (positions[u] < positions[v] && x > 0.5) {
        return u_before_v;
    } else if (positions[v] < positions[u] & x <= 0.5) {
        return v_before_u;
    } else {
        return indeterminate;
    }
}

pattern instance::based_on_pattern(const vertex_t &u,              //
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

};  // namespace pace
