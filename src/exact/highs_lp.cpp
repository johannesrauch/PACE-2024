#include "exact/highs_lp.hpp"

namespace pace {

//
// ctor
//

highs_lp::highs_lp(instance &instance_, highs_lp_params params)
    : highs_base(instance_),
      info{u_old, v_old, w_old},
      params(params),
      sort_h(instance_) {
    // this->params.max_new_rows =
    //     std::max(static_cast<std::size_t>(0.2 * std::sqrt(n_free) * n_free),
    //              params.max_new_rows);

    rows_to_delete.reserve(params.max_new_rows);
    lower_bounds.reserve(params.max_new_rows);
    upper_bounds.reserve(params.max_new_rows);
    starts.reserve(params.max_new_rows);
    indices.reserve(3 * params.max_new_rows);
    values.reserve(3 * params.max_new_rows);

    assert(n_free >= 2);
    this->params.max_new_rows = n_free_2 * (n_free - 2) / (3 * 50);

    info.n_cols = get_n_cols();
    // PACE_DEBUG_PRINTF("start add_initial_rows\n");
    // add_some_interesting_rows_based_on_ordering();
    // if (get_n_rows() == 0) add_some_interesting_rows();
    // PACE_DEBUG_PRINTF("end   add_initial_rows\n");
}

//
// public methods
//

//
// lp solving methods
//

bool highs_lp::cut() {
    n_old_rows = get_n_rows();
    bool success = check_3cycles_backarcs();
    return success;
}

void highs_lp::run(const int32_t max_simplex_iter, const bool bookkeeping) {
    set_simplex_iteration_limit(max_simplex_iter);
    solver.run();
    if (bookkeeping) update_simplex_info();
}

void highs_lp::initial_solve() {
    info.n_solve_iters = 0;
    const std::size_t max_new_rows = params.max_new_rows;
    // params.max_new_rows = std::numeric_limits<std::size_t>::max();
    // add_all_interesting_rows();
    PACE_DEBUG_PRINTF_LPINFO_LINE();

    for (;;) {
        run();
        PACE_DEBUG_PRINTF_LPINFO(info);
        ++info.n_solve_iters;
        reset_row_info();
        if (get_rounded_objective_value() >= upper_bound) break;
        if (!cut()) break;
        delete_positive_slack_rows();
    }

    params.max_new_rows = max_new_rows;
}

void highs_lp::initial_partial_solve() {
    info.n_solve_iters = 0;
    info.n_3cycle_iters = 0;
    PACE_DEBUG_PRINTF_LPINFO_LINE();

    // few 3-cycle iterations, solved to optimality, to sieve for "important"
    // 3-cycle ieqs
    bool cut_generated{true};
    while (cut_generated &&
           info.n_3cycle_iters < params.max_initial_solve_3cycle_iters &&
           info.n_solve_iters < params.max_initial_solve_iters) {
        if (info.was_warmstart &&
            info.n_simplex_iters > info.n_simplex_coldstart_iters / 2) {
            delete_positive_slack_rows();
        }
        run(params.max_initial_solve_simplex_iters);

        PACE_DEBUG_PRINTF_LPINFO(info);
        ++info.n_solve_iters;
        reset_row_info();

        if (get_rounded_objective_value() >= upper_bound) return;
        cut_generated = cut();
    }
}

void highs_lp::resolve() {
    info.n_solve_iters = 0;
    const std::size_t max_new_rows = params.max_new_rows;
    params.max_new_rows = std::numeric_limits<std::size_t>::max();
    PACE_DEBUG_PRINTF_LPINFO_LINE();

    do {
        run();
        PACE_DEBUG_PRINTF_LPINFO(info);
        ++info.n_solve_iters;
        reset_row_info();
        if (get_rounded_objective_value() >= upper_bound) break;
    } while (cut());

    params.max_new_rows = max_new_rows;
}

//
// row modification methods
//

void highs_lp::delete_positive_slack_rows() {
    info.tried_deleting_rows = true;

    // gather rows to delete
    rows_to_delete.clear();
    auto it_rows = rows.begin();
    for (std::size_t i = 0; i < n_old_rows; ++i) {
        assert(it_rows != rows.end());
        auto it = rows_info.find(*it_rows);
        assert(it != rows_info.end());

        const bool has_slack = has_row_slack(i);
        const bool spare = it->second.n_deleted > it->second.n_spared;
        const bool remove = has_slack && !spare;

        if (remove) {
            rows_to_delete.emplace_back(i);
            ++it->second.n_deleted;
            it->second.n_spared = 0;
            it_rows = rows.erase(it_rows);
        } else {
            ++it_rows;
        }

        if (has_slack && spare) {
            ++it->second.n_spared;
            ++info.n_delete_rows_spared;
        }
    }

    info.n_deleted_rows = rows_to_delete.size();
    if (info.n_deleted_rows > 0) {
        solver.deleteRows(info.n_deleted_rows, &rows_to_delete[0]);
    }
}

//
// column modification methods
//

void highs_lp::fix_columns() {
    for (std::size_t j = 0; j < get_n_cols(); ++j) {
        const crossing_number_t diff = upper_bound - lower_bound();
        const double coeff = get_objective_coefficient(j);

        if (std::abs(coeff) >= static_cast<double>(diff) && coeff != 0.) {
            fix_column(j, coeff > 0 ? 0. : 1.);
            PACE_DEBUG_PRINTF("(optimality condition)\n");
        }
    }
}

//
// information methods
//

bool highs_lp::is_integral() {
    for (std::size_t j = 0; j < get_n_cols(); ++j) {
        if (!is_column_integral(j)) {
            return false;
        }
    }
    return true;
}

bool highs_lp::is_3cycle_interesting(const vertex_t &u, const vertex_t &v,
                                     const vertex_t &w) {
    assert(u < v);
    assert(v < w);

    const crossing_matrix &cr = cr_matrix();

    const crossing_number_t c_uvw = cr(u, v) + cr(u, w) + cr(v, w);
    const crossing_number_t c_uwv = cr(u, v) + cr(u, w) + cr(w, v);
    const crossing_number_t c_vuw = cr(v, u) + cr(u, w) + cr(v, w);
    const crossing_number_t c_vwu = cr(v, u) + cr(w, u) + cr(v, w);
    const crossing_number_t c_wuv = cr(u, v) + cr(w, u) + cr(w, v);
    const crossing_number_t c_wvu = cr(v, u) + cr(w, u) + cr(w, v);
    const crossing_number_t c_min_allowed =
        std::min({c_uvw, c_uwv, c_vuw, c_vwu, c_wuv, c_wvu});

    const crossing_number_t c_uv_vw_wu = cr(u, v) + cr(v, w) + cr(w, u);
    const crossing_number_t c_vu_wv_uw = cr(v, u) + cr(w, v) + cr(u, w);
    const crossing_number_t c_min_forbidden = std::min(c_uv_vw_wu, c_vu_wv_uw);

    return c_min_forbidden < c_min_allowed;
}

//
// 3-cycle methods
//

bool highs_lp::add_3cycle_row_to_internal_rows(  //
    const vertex_t &u, const vertex_t &v, const vertex_t &w) {
    assert(u < v);
    assert(v < w);
    rows.emplace_back(u, v, w);
    const triple uvw(u, v, w);
    auto it = rows_info.find(uvw);
    if (it == rows_info.end()) {
        rows_info.emplace(uvw, row_info());
        return true;
    } else {
        return false;
    }
}

void highs_lp::add_some_interesting_rows_based_on_ordering() {
    clear_aux_vectors();
    const std::vector<vertex_t> &ordering = get_ordering();
    std::vector<vertex_t> positions(ordering.size());
    inverse(ordering, positions);
    const crossing_matrix &cr = cr_matrix();

    std::vector<std::pair<vertex_t, vertex_t>> unnatural;
    for (const auto &[u, v] : unsettled_pairs()) {
        const crossing_number_t &c_uv = cr(u, v);
        const crossing_number_t &c_vu = cr(v, u);
        if (c_uv < c_vu && positions[u] > positions[v] ||
            c_uv > c_vu && positions[u] < positions[v]) {
            unnatural.emplace_back(u, v);
        }
    }

    const std::size_t n_unnatural = unnatural.size();
    if (n_unnatural == 0) return;
    const std::size_t n_rows_pre_unnatural =
        (params.max_initial_rows + n_unnatural - 1) / n_unnatural;
    std::uniform_int_distribution<vertex_t> distribution(0, n_free - 1);

    for (const auto &[u, v] : unnatural) {
        std::unordered_set<vertex_t> set;
        for (std::size_t i = 0; i < n_rows_pre_unnatural; ++i) {
            set.insert(distribution(rd_generator));
        }
        set.erase(u);
        set.erase(v);
        assert(u < v);
        for (const vertex_t &w : set) {
            if (w < u) {
                assert(w < v);
                if (get_n_vars_in_lp(w, u, v) >= 2 &&
                    add_3cycle_row_to_internal_rows(w, u, v))
                    add_3cycle_row_to_aux_vectors(w, u, v);
            } else if (w < v) {
                assert(u < w);
                if (get_n_vars_in_lp(u, w, v) >= 2 &&
                    add_3cycle_row_to_internal_rows(u, w, v))
                    add_3cycle_row_to_aux_vectors(u, w, v);
            } else {
                assert(v < w);
                if (get_n_vars_in_lp(u, v, w) >= 2 &&
                    add_3cycle_row_to_internal_rows(u, v, w))
                    add_3cycle_row_to_aux_vectors(u, v, w);
            }
        }
    }

    solver.addRows(lower_bounds.size(), &lower_bounds[0], &upper_bounds[0],  //
                   indices.size(), &starts[0], &indices[0], &values[0]);
    info.n_rows = get_n_rows();
}

void highs_lp::add_some_interesting_rows() {
    clear_aux_vectors();
    std::vector<triple> candidates;
    for (vertex_t u = 0u; u + 2u < n_free; ++u) {
        for (vertex_t v = u + 1u; v + 1u < n_free; ++v) {
            for (vertex_t w = v + 1u; w < n_free; ++w) {
                assert(u < v);
                assert(v < w);
                if (get_n_vars_in_lp(u, v, w) >= 2 &&
                    is_3cycle_interesting(u, v, w)) {
                    candidates.emplace_back(u, v, w);
                }
            }
        }
    }

    const double p = static_cast<double>(params.max_initial_rows) -
                     static_cast<double>(get_n_rows()) / candidates.size();
    for (const auto &[u, v, w] : candidates) {
        if (coinflip(p)) {
            add_3cycle_row_to_aux_vectors(u, v, w);
            add_3cycle_row_to_internal_rows(u, v, w);
        }
    }

    solver.addRows(lower_bounds.size(), &lower_bounds[0],
                   &upper_bounds[0],  //
                   indices.size(), &starts[0], &indices[0], &values[0]);
    info.n_init_rows_candidates = candidates.size();
    info.n_rows = get_n_rows();
}

void highs_lp::add_all_interesting_rows() {
    const std::vector<std::pair<vertex_t, vertex_t>> unsettled =
        unsettled_pairs();

    auto end = unsettled.end();
    for (auto it1 = unsettled.begin(); it1 + 1 != end; ++it1) {
        const auto &[u, v] = *it1;
        for (auto it2 = it1 + 1; it2 != end; ++it2) {
            const auto &[x, y] = *it2;
            if (u == x) {
                if (is_3cycle_interesting(u, v, y)) {
                    add_3cycle_row_to_aux_vectors(u, v, y);
                    add_3cycle_row_to_internal_rows(u, v, y);
                }
            } else
                break;
        }

        for (auto it2 = std::lower_bound(it1, end, std::make_pair(v, 0u));
             it2 != end; ++it2) {
            const auto &[x, y] = *it2;
            if (v == x) {
                if (is_3cycle_interesting(u, v, y)) {
                    add_3cycle_row_to_aux_vectors(u, v, y);
                    add_3cycle_row_to_internal_rows(u, v, y);
                }
            } else
                break;
        }
    }
}

std::size_t highs_lp::add_3cycle_rows() {
    info.n_added_rows = std::min(get_n_bucket_entries(),
                                 static_cast<std::size_t>(params.max_new_rows));
    if (info.n_added_rows <= 0) return 0;

    clear_aux_vectors();

    std::size_t i = 0;
    bool go_on = true;
    for (auto r_it = buckets.rbegin(); r_it != buckets.rend() && go_on;
         ++r_it) {
        for (const auto &[u, v, w] : *r_it) {
            add_3cycle_row_to_aux_vectors(u, v, w);
            add_3cycle_row_to_internal_rows(u, v, w);
            ++i;
            if (i >= params.max_new_rows) {
                go_on = false;
                break;
            }
        }
    }
    assert(i == info.n_added_rows);

    solver.addRows(info.n_added_rows, &lower_bounds[0],
                   &upper_bounds[0],  //
                   indices.size(), &starts[0], &indices[0], &values[0]);
    return info.n_added_rows;
}

bool highs_lp::check_3cycle(const vertex_t &u, const vertex_t &v,
                            const vertex_t &w) {
    assert(u < v);
    assert(v < w);
    const auto [uv, vw, uw] = flat_indices(n_free, n_free_2, u, v, w);
    const double x = get_3cycle_value(uv, vw, uw);
    const double interval_width = 1. + 3. * params.tol_feasibility;
    if (is_3cycle_lb_violated(x)) {
        const double x_normalized = -x / interval_width;
        get_bucket(x_normalized).emplace_back(u, v, w);
        // get_bucket_from_score(u, v, w).emplace_back(u, v, w);
        return true;
    }
    if (is_3cycle_ub_violated(x)) {
        const double ub = 1. + params.tol_feasibility;
        const double x_normalized = (x - ub) / interval_width;
        get_bucket(x_normalized).emplace_back(u, v, w);
        // get_bucket_from_score(u, v, w).emplace_back(u, v, w);
        return true;
    }
    return false;
}

bool highs_lp::check_3cycles() {
    clear_buckets();
    bool stop = false;
    info.new_3cycle_iter = false;
    // yes, I use the dark forces here, because it speeds things up and is
    // imo cleaner
    vertex_t u = u_old, v = v_old, w = w_old;
    goto check_3cycles_in_for;
    for (; u + 2u < n_free; ++u) {
        for (v = u + 1; v + 1u < n_free; ++v) {
            for (w = v + 1; w < n_free; ++w) {
            check_3cycles_in_for:  // to pick off where we left
                check_3cycle(u, v, w);
                stop = is_last_bucket_full() || is_n_bucket_entries_large();
                if (stop) goto check_3cycles_after_for;
            }
        }
    }
    update_3cycle_iteration_info();
    for (u = 0; u + 2u < n_free; ++u) {
        bool u_eq_old = u == u_old;
        for (v = u + 1; v + 1u < n_free; ++v) {
            bool v_eq_old = v == v_old;
            for (w = v + 1; w < n_free; ++w) {
                if (u_eq_old && v_eq_old && w == w_old)
                    goto check_3cycles_after_for;
                check_3cycle(u, v, w);
                stop = is_last_bucket_full() || is_n_bucket_entries_large();
                if (stop) goto check_3cycles_after_for;
            }
        }
    }

check_3cycles_after_for:
    info.n_bucket_entries = get_n_bucket_entries();
    u_old = u, v_old = v, w_old = w;
    return add_3cycle_rows() > 0;
}

bool highs_lp::check_3cycles_depr() {
    clear_buckets();
    const std::vector<std::pair<vertex_t, vertex_t>> unsettled =
        unsettled_pairs();

    auto end = unsettled.end();
    bool stop = false;
    for (auto it1 = unsettled.begin(); it1 + 1 != end; ++it1) {
        const auto &[u, v] = *it1;
        for (auto it2 = it1 + 1; it2 != end; ++it2) {
            const auto &[x, y] = *it2;
            if (u == x) {
                check_3cycle(u, v, y);
                stop = is_last_bucket_full() || is_n_bucket_entries_large();
                if (stop) goto check_3cycles_fast_after_for;
            } else
                break;
        }

        for (auto it2 = std::lower_bound(it1, end, std::make_pair(v, 0u));
             it2 != end; ++it2) {
            const auto &[x, y] = *it2;
            if (v == x) {
                check_3cycle(u, v, y);
                stop = is_last_bucket_full() || is_n_bucket_entries_large();
                if (stop) goto check_3cycles_fast_after_for;
            } else
                break;
        }
    }

check_3cycles_fast_after_for:
    info.n_bucket_entries = get_n_bucket_entries();
    return add_3cycle_rows() > 0;
}

bool highs_lp::check_3cycles_backarcs() {
    clear_buckets();

    std::vector<vertex_t> ordering;
    sort_h(*this, ordering);
    std::vector<vertex_t> positions;
    inverse(ordering, positions);

    for (std::size_t i = 0; i + 2u < n_free; ++i) {
        const vertex_t a = positions[i];
        for (std::size_t j = i + 1u; j + 1u < n_free; ++j) {
            const vertex_t b = positions[j];
            assert(a != b);

            vertex_t u, v;
            bool backarc;
            if (a < b) {
                u = a, v = b;
                backarc = get_variable_value(u, v) < 1. - params.tol_integer;
            } else {
                u = b, v = a;
                backarc = get_variable_value(u, v) > params.tol_integer;
            }

            if (backarc) {
                for (std::size_t k = i + 1u; k < j; ++k) {
                    const vertex_t w = positions[k];
                    if (a < w &&
                        get_variable_value(a, w) < 1. - params.tol_integer)
                        continue;
                    if (w < a && get_variable_value(w, a) > params.tol_integer)
                        continue;
                    
                    if (w < u)
                        check_3cycle(w, u, v);
                    else if (w < v)
                        check_3cycle(u, w, v);
                    else
                        check_3cycle(u, v, w);
                }

                for (std::size_t k = j + 1u; k < n_free; ++k) {
                    const vertex_t w = positions[k];
                    if (w < u)
                        check_3cycle(w, u, v);
                    else if (w < v)
                        check_3cycle(u, w, v);
                    else
                        check_3cycle(u, v, w);
                }
            }
        }
    }

    info.n_bucket_entries = get_n_bucket_entries();
    return add_3cycle_rows() > 0;
}

//
// private methods
//

//
// bookkeeping methods
//

void highs_lp::update_3cycle_iteration_info() {
    if (info.n_3cycle_iters < params.max_initial_solve_3cycle_iters) {
        params.max_new_rows *= 2;
    }
    ++info.n_3cycle_iters;
    info.new_3cycle_iter = true;
}

void highs_lp::update_simplex_info() {
    info.n_simplex_iters = solver.getSimplexIterationCount();
    if (info.n_simplex_coldstart_iters == 0 || info.n_deleted_rows > 0) {
        info.n_simplex_coldstart_iters = info.n_simplex_iters;
    }
    info.was_warmstart = info.n_deleted_rows == 0;
    info.n_avg_simplex_iters =
        (info.n_avg_simplex_iters * info.n_runs + info.n_simplex_iters) /
        (info.n_runs + 1);

    info.n_rows = get_n_rows();
    ++info.n_runs;

    info.t_simplex = solver.getRunTime();
    info.objective_value = get_objective_value();
    info.percent_integral = 0.;
    for (std::size_t j = 0; j < get_n_cols(); ++j)
        info.percent_integral += is_column_integral(j);
    info.percent_integral /= get_n_cols();
}

void highs_lp::reset_row_info() {
    info.n_bucket_entries = 0;
    info.n_added_rows = 0;
    info.n_deleted_rows = 0;
    info.n_delete_rows_spared = 0;
    info.tried_deleting_rows = false;
}

}  // namespace pace
