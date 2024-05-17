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

instance *instance::new_rins_instance(highs_lp &lp) {
    instance *subinstance = new_subinstance();
    assert(!subinstance->restriction_graph_ptr);
    subinstance->create_kernel(&lp);
    return subinstance;
}

instance *instance::new_lsearch_instance(const uint16_t lsearch_width) {
    assert(lsearch_width > 0);
    instance *subinstance = new_subinstance();
    assert(!subinstance->restriction_graph_ptr);
    subinstance->create_kernel(nullptr, lsearch_width);
    return subinstance;
}

instance *instance::new_subinstance() {
    instance *subinstance = new instance(graph);
    subinstance->cr_matrix_ptr =
        const_cast<crossing_matrix *>(&get_cr_matrix());
    subinstance->lower_bound = get_lower_bound();
    subinstance->upper_bound = get_upper_bound();
    subinstance->ordering = get_ordering();
    return subinstance;
}

//
// private methods
//

void pace::instance::create_kernel(highs_lp *lp, const uint16_t lsearch_width) {
    if (restriction_graph_ptr) return;

    PACE_DEBUG_PRINTF("start create_kernel\n");
    get_cr_matrix();  // initializes crossing matrix
    const std::size_t n_free = graph.get_n_free();
    restriction_graph_ptr = std::make_unique<digraph>(graph.get_n_free());
    magic.resize(n_free_2);
    std::vector<vertex_t> positions;
    if (lp != nullptr || lsearch_width > 0) inverse(ordering, positions);

    // let's foresee if we can fix u < v or v < u
    for (vertex_t u = 0u; u + 1u < n_free; ++u) {
        for (vertex_t v = u + 1u; v < n_free; ++v) {
            switch (foresee(u, v, lp, positions, lsearch_width)) {
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
                          const std::vector<vertex_t> &positions,
                          const uint16_t &lsearch_width) {
    assert(u < v);

    pattern p = indeterminate;
    if (lp != nullptr) {
        p = based_on_relaxation(u, v, *lp, positions);
        if (p != indeterminate) return p;
    }

    if (lsearch_width > 0) {
        p = based_on_ordering(u, v, positions, lsearch_width);
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

    if (positions[u] < positions[v] && x > 0.5) return u_before_v;
    if (positions[v] < positions[u] & x <= 0.5) return v_before_u;
    return indeterminate;
}

pattern instance::based_on_ordering(const vertex_t &u, const vertex_t &v,
                                    const std::vector<vertex_t> &positions,
                                    const uint16_t &lsearch_width) {
    if (positions[u] + lsearch_width < positions[v]) return u_before_v;
    if (positions[v] + lsearch_width < positions[u]) return v_before_u;
    return indeterminate;
}

pattern instance::based_on_degree(const vertex_t &u, const vertex_t &v) {
    if (graph.get_degree(u) == 0) return pattern::u_before_v;
    if (graph.get_degree(v) == 0) return pattern::v_before_u;
    return pattern::indeterminate;
}

pattern instance::based_on_crossing_numbers(const crossing_number_t &c_uv,
                                            const crossing_number_t &c_vu) {
    assert(c_uv != c_vu);
    if (c_uv == 0) return pattern::u_before_v;
    if (c_vu == 0) return pattern::v_before_u;
    return pattern::indeterminate;
}

pattern instance::based_on_bounds(const crossing_number_t &c_uv,
                                  const crossing_number_t &c_vu) {
    const crossing_number_t diff = upper_bound - lower_bound;
    if (c_uv > c_vu && c_uv - c_vu > diff) return pattern::v_before_u;
    if (c_vu > c_uv && c_vu - c_uv > diff) return pattern::u_before_v;
    return pattern::indeterminate;
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
