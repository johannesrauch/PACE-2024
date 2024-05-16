#include "exact/branch_and_cut.hpp"

namespace pace {

//
// private methods
//

void branch_and_cut::build_ordering(std::vector<vertex_t> &ordering) {
    assert(lp_solver_ptr->is_integral());
#ifndef NDEBUG
    const bool acyclic =
#endif
        build_restr_graph_ordering(  //
            lp_solver_ptr->get_column_values(), unsettled_pairs(),
            restriction_graph(), ordering);
    assert(acyclic);

    crossing_number_t n_crossings =
        lp_solver_ptr->get_rounded_objective_value();
    assert(n_crossings < upper_bound);
    assert(n_crossings == number_of_crossings(graph, ordering));

    // try to improve new solution
    n_crossings = shift_h(ordering, n_crossings);
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
            lp_solver_ptr->fix_column(
                node.j,
                lp_solver_ptr->get_column_value(node.j) > 0.5 ? 0. : 1.);
            PACE_DEBUG_PRINTF("(backtrack)\n");
            node.fix_opposite = false;
            return false;
        } else {
            lp_solver_ptr->unfix_column(node.j);
            PACE_DEBUG_PRINTF("(backtrack)\n");
            stack.pop();
        }
    }
    return true;
}

void branch_and_cut::branch() {
    info.n_iter_3cycle_current_node = 0;

    // gather branch node info
    const double obj_val_old = lp_solver_ptr->get_objective_value();
    const std::size_t j = reli_branch_ptr->get_branching_column();
    const double x_j_old = reli_branch_ptr->get_column_value(j);

    // branch
    ++info.n_branch_nodes;
    lp_solver_ptr->fix_column(
        j, reli_branch_ptr->get_up_cost(j) > reli_branch_ptr->get_down_cost(j));
    PACE_DEBUG_PRINTF("(branch)\n");
    stack.emplace(branch_node{j, obj_val_old, x_j_old, true});
}

bool branch_and_cut::branch_and_bound_and_cut(std::vector<vertex_t> &ordering) {
    lp_solver_ptr->resolve();

    if (lp_solver_ptr->is_infeasible()) {
        return backtrack();
    }

    if (lp_solver_ptr->get_rounded_objective_value() >= upper_bound) {
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
    if (heuristic.informed(*lp_solver_ptr, another_ordering) < ub_old)
        std::swap(ordering, another_ordering);
    info.relax_h_confidence = heuristic.get_confidence();

    branch();
    return false;
}

//
// public methods
//

uint32_t branch_and_cut::operator()(std::vector<vertex_t> &ordering) {
    PACE_DEBUG_PRINTF("\nstart branch and cut\n");

    // heuristics
    info.n_crossings_h = heuristic.uninformed(ordering);
    if (lower_bound() >= upper_bound) return upper_bound;

    // lp
    if (!lp_solver_ptr) lp_solver_ptr = std::make_unique<highs_lp>(instance_);
    if (unsettled_pairs().size() == 0) {
        build_ordering(ordering);
        return upper_bound;
    }
    if (!reli_branch_ptr)
        reli_branch_ptr =
            std::make_unique<reliability_branching>(*lp_solver_ptr);

    const highs_lp_info &info_lp = lp_solver_ptr->get_info();
    PACE_DEBUG_PRINTF("%11s=%11u, %11s=%11u, %11s=%11u\n",  //
                      "n cols", info_lp.n_cols,             //
                      "n initrows", info_lp.n_rows,         //
                      "n cand", info_lp.n_init_rows_candidates);

    // driver
    lp_solver_ptr->initial_partial_solve();
    while (!branch_and_bound_and_cut(ordering)) {
        ++info.n_iters;
    };

    PACE_DEBUG_PRINTF("end   branch and cut\n");
    return upper_bound;
}

};  // namespace pace
