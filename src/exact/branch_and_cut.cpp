#include "exact/branch_and_cut.hpp"

namespace pace {

//
// private methods
//

void branch_and_cut::build_ordering(std::vector<vertex_t> &ordering) {
    assert(lp_solver_ptr->is_integral());
#ifndef NDEBUG
    std::vector<std::pair<vertex_t, vertex_t>> backarcs;
    digraph &restr_graph = restriction_graph();
    build_restr_graph_ordering(  //
        lp_solver_ptr->get_column_values(), unsettled_pairs(),
        lp_solver_ptr->get_tol_integer(), restr_graph,  //
        ordering, backarcs);
    if (backarcs.size() > 0) {
        PACE_DEBUG_PRINTF("ordering:\n");
        PACE_DEBUG_PRINTF_VECTOR(ordering);
        PACE_DEBUG_PRINTF("backarcs:\n");
        PACE_DEBUG_PRINTF_VECTOR(backarcs);
        const auto &[a, b] = backarcs[0];
        for (const vertex_t u : restr_graph.get_neighbors(b)) {
            for (const vertex_t v : restr_graph.get_neighbors(u)) {
                if (v == a) {
                    PACE_DEBUG_PRINTF("%u->%u->%u->%u\n", a, b, u, v);
                    std::vector<vertex_t> abu = {a, b, u};
                    std::sort(abu.begin(), abu.end());
                    const auto [xy, yz, xz] =
                        flat_indices(n_free, n_free_2, abu[0], abu[1], abu[2]);
                    PACE_DEBUG_PRINTF("xy=%1.1f, yz=%1.1f, xz=%1.1f\n",
                                      lp_solver_ptr->get_variable_value(xy),
                                      lp_solver_ptr->get_variable_value(yz),
                                      lp_solver_ptr->get_variable_value(xz));
                }
            }
        }
    }
    assert(backarcs.size() == 0);
#else
    build_restr_graph_ordering(
        lp_solver_ptr->get_column_values(), unsettled_pairs(),
        lp_solver_ptr->get_tol_integer(), restriction_graph(), ordering);
#endif

    crossing_number_t n_crossings =
        lp_solver_ptr->get_rounded_objective_value();
    assert(n_crossings < upper_bound);
    assert(n_crossings == number_of_crossings(graph, ordering));

    // try to improve new solution
    n_crossings = heuristic.shift(ordering, n_crossings);
    update_ordering(ordering, n_crossings);
}

void branch_and_cut::update_costs() {
    if (stack.empty()) return;
    const auto &[j, obj_val_old, x_j_old, _] = stack.top();
    const double obj_val_new = lp_solver_ptr->get_objective_value();
    if (lp_solver_ptr->get_column_value(j) > 0.5) {
        reli_branch_ptr->update_up_cost(j, obj_val_new, obj_val_old, x_j_old);
    } else {
        reli_branch_ptr->update_down_cost(j, obj_val_new, obj_val_old, x_j_old);
    }
}

bool branch_and_cut::backtrack() {
    info.n_iter_3cycle_current_node = 0;

    while (!stack.empty()) {
        branch_node &node = stack.top();
        if (node.fix_opposite) {
            ++info.n_search_nodes;
            lp_solver_ptr->fix_column(
                node.j,
                lp_solver_ptr->get_column_value(node.j) > 0.5 ? 0. : 1.);
            PACE_DEBUG_PRINTF("(backtrack)\n");
            node.fix_opposite = false;
            return info.n_search_nodes > params.max_nodes;
        } else {
            lp_solver_ptr->unfix_column(node.j);
            PACE_DEBUG_PRINTF("(backtrack)\n");
            stack.pop();
        }
    }

    return true;
}

bool branch_and_cut::branch() {
    info.n_iter_3cycle_current_node = 0;

    // gather branch node info
    const double obj_val_old = lp_solver_ptr->get_objective_value();
    const std::size_t j = reli_branch_ptr->get_branching_column();
    const double x_j_old = reli_branch_ptr->get_column_value(j);

    // branch
    ++info.n_search_nodes;
    lp_solver_ptr->fix_column(
        j, reli_branch_ptr->get_up_cost(j) > reli_branch_ptr->get_down_cost(j));
    PACE_DEBUG_PRINTF("(branch)\n");
    stack.emplace(branch_node{j, obj_val_old, x_j_old, true});

    return info.n_search_nodes > params.max_nodes;
}

bool branch_and_cut::branch_and_bound_and_cut(std::vector<vertex_t> &ordering) {
    if (stack.empty())
        lp_solver_ptr->initial_solve();
    else
        lp_solver_ptr->resolve();

    if (lp_solver_ptr->is_infeasible()) {
        return backtrack();
    }

    const crossing_number_t rounded_obj_val =
        lp_solver_ptr->get_rounded_objective_value();
    if (rounded_obj_val >= upper_bound) {
        update_costs();
        return backtrack();
    }

    if (lp_solver_ptr->is_integral()) {
        build_ordering(ordering);
        lp_solver_ptr->fix_columns();
        update_costs();
        return backtrack();
    }

    const crossing_number_t ub_old = upper_bound;
    const bool do_rins =
        params.do_rins && info.n_search_nodes % params.n_nodes_until_rins == 0;
    heuristic.informed(*lp_solver_ptr, ordering, do_rins);
    info.relax_h_confidence = heuristic.get_confidence();

    if (rounded_obj_val >= upper_bound) {
        update_costs();
        return backtrack();
    }

    if (upper_bound < ub_old) lp_solver_ptr->fix_columns();

    return branch();
}

//
// public methods
//

uint32_t branch_and_cut::operator()(std::vector<vertex_t> &ordering) {
    PACE_DEBUG_PRINTF("start branch and cut\n");
    info.t_start = now();
    cr_matrix();  // initialize

    // uninformed heuristics
    if (params.do_uninformed_h) {
        info.n_crossings_h = heuristic.uninformed(ordering, params.do_lsearch);
    }
    PACE_DEBUG_PRINTF_BOUNDS(lower_bound(), upper_bound);
    if (lower_bound() >= upper_bound) return upper_bound;

    // create lp
    if (!lp_solver_ptr) lp_solver_ptr = std::make_unique<highs_lp>(instance_);
    if (unsettled_pairs().size() == 0) {
        build_ordering(ordering);
        return upper_bound;
    }
    if (!reli_branch_ptr)
        reli_branch_ptr =
            std::make_unique<reliability_branching>(*lp_solver_ptr);

    // uninformed rins heuristics
    if (params.do_uninformed_rins) {
        lp_solver_ptr->run();
        heuristic.rins(*lp_solver_ptr, ordering);
    }
    if (lower_bound() >= upper_bound) return upper_bound;

    const highs_lp_info &info_lp = lp_solver_ptr->get_info();
    PACE_DEBUG_PRINTF("%11s=%11u, %11s=%11u, %11s=%11u\n",  //
                      "n cols", info_lp.n_cols,             //
                      "n initrows", info_lp.n_rows,         //
                      "n cand", info_lp.n_init_rows_candidates);

    // driver
    if (params.do_initial_partial_solve) {
        PACE_DEBUG_PRINTF("start initial partial solve\n");
        lp_solver_ptr->initial_partial_solve();
        PACE_DEBUG_PRINTF("end   initial partial solve\n");
    }
    while (!branch_and_bound_and_cut(ordering)) {
        ++info.n_iters;
        PACE_DEBUG_PRINTF_BOUNDS(lower_bound(), upper_bound);
    };

    PACE_DEBUG_PRINTF("end   branch and cut\n");
    info.n_rows = lp_solver_ptr->get_n_rows();
    info.t_end = now();
    info.t_sol = instance_.get_t_sol();
    PACE_DEBUG_PRINTF_BOUNDS(lower_bound(), upper_bound);
    PACE_DEBUG_PRINTF_SUMMARY(info);

    ordering = get_ordering();
    return upper_bound;
}

};  // namespace pace
