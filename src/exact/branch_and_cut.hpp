#ifndef PACE_EXACT_BRANCH_AND_CUT_HPP
#define PACE_EXACT_BRANCH_AND_CUT_HPP

#include <functional>
#include <stack>
#include <unordered_map>

#include "exact/highs_lp.hpp"
#include "exact/reliability_branching.hpp"
#include "heuristic/heuristics.hpp"
#include "io/output.hpp"
#include "log/debug_printf.hpp"
#include "model/instance.hpp"
#include "utils/crossings_utils.hpp"
#include "utils/restr_graph_utils.hpp"

namespace pace {

struct branch_and_cut_params {
    bool delete_rows{true};
    std::size_t max_iter_3cycle_until_branch{1000};
};

class branch_and_cut : public instance_view {
    /**
     * @brief relevant information of pre-branching situation
     */
    struct branch_node {
        const std::size_t j;
        const double obj_val_old;
        const double x_j_old;
        bool fix_opposite;
        // HighsInt frozen_basis_id;
    };

    /**
     * @brief for a dfs through search tree
     */
    std::stack<branch_node> stack;

    /**
     * @brief interface to the ilp relaxation solver
     */
    std::unique_ptr<highs_lp> lp_solver_ptr;
    std::unique_ptr<reliability_branching> reli_branch_ptr;

    branch_and_cut_info info;
    branch_and_cut_params params;

    shift_heuristic shift_h;
    heuristics heuristic;
    std::vector<vertex_t> another_ordering;

   public:
    /**
     * @brief constructs and initializes the branch and cut solver
     */
    branch_and_cut(instance &instance_, branch_and_cut_params params = branch_and_cut_params())
        : instance_view(instance_),
          info{lower_bound(), upper_bound},
          params(params),
          shift_h(instance_),
          heuristic(instance_),
          another_ordering(n_free) {
        assert(n_free > 0);
    }

    // delete copy constructor, move constructor, copy assignment and move assignment
    branch_and_cut(const branch_and_cut &rhs) = delete;
    branch_and_cut(branch_and_cut &&rhs) = delete;
    branch_and_cut &operator=(const branch_and_cut &rhs) = delete;
    branch_and_cut &operator=(branch_and_cut &&rhs) = delete;

   private:
    /**
     * @brief constructs the restrictions graph into restriction_graph.
     * the topological sort of `restriction_graph` gives an ordering of the vertices of the free
     * layer of `graph`, which we store in `ordering`.
     */
    void build_ordering(std::vector<vertex_t> &ordering) {
        assert(lp_solver_ptr->is_integral());
#ifndef NDEBUG
        const bool acyclic =
#endif
            build_restr_graph_ordering(  //
                lp_solver_ptr->get_column_values(), unsettled_pairs(), restriction_graph(), ordering);
        assert(acyclic);

        crossing_number_t n_crossings = lp_solver_ptr->get_rounded_objective_value();
        assert(n_crossings < upper_bound);
        assert(n_crossings == number_of_crossings(graph, ordering));

        // try to improve new solution
        n_crossings = shift_h(ordering, n_crossings);
        update_ordering(ordering, n_crossings);
    }

    /**
     * @brief why u not work?!
     */
    void _build_ordering(std::vector<vertex_t> ordering) {
        assert(lp_solver_ptr->is_integral());
        ordering.resize(n_free);
        for (vertex_t u = 0; u < n_free; ++u) {
            std::size_t i = 0;
            for (vertex_t v = 0; v < u; ++v) {
                i += lp_solver_ptr->get_variable_value(flat_index(n_free, n_free_2, v, u)) > 0.5;
            }
            for (vertex_t v = u + 1u; v < n_free; ++v) {
                i += lp_solver_ptr->get_variable_value(flat_index(n_free, n_free_2, u, v)) < 0.5;
            }
            assert(i < n_free);
            ordering[i] = u;
        }
        assert(test::is_permutation(ordering));

        crossing_number_t n_crossings = lp_solver_ptr->get_rounded_objective_value();
        assert(n_crossings < upper_bound);
        assert(n_crossings == number_of_crossings(graph, ordering));
        n_crossings = shift_h(ordering, n_crossings);
        update_ordering(ordering, n_crossings);
    }

    /**
     * @brief updates pseudo-costs
     */
    void update_costs() {
        if (stack.empty()) return;
        const auto &[j, obj_val_old, x_j_old, _] = stack.top();
        const double obj_val_new = lp_solver_ptr->get_objective_value();
        if (lp_solver_ptr->get_column_value(j) > 0.5) {
            reli_branch_ptr->update_up_cost(j, obj_val_new, obj_val_old, x_j_old);
        } else {
            reli_branch_ptr->update_down_cost(j, obj_val_new, obj_val_old, x_j_old);
        }
    }

    //
    // branch and bound and cut methods
    //

    /**
     * @brief if there's a variable on the stack,
     * the function pops it, and we fix the opposite value or we backtrack again.
     *
     * @return true if stack is empty
     * @return false otherwise
     */
    bool backtrack() {
        info.n_iter_3cycle_current_node = 0;
        while (!stack.empty()) {
            branch_node &node = stack.top();
            if (node.fix_opposite) {
                lp_solver_ptr->fix_column(node.j, lp_solver_ptr->get_column_value(node.j) > 0.5 ? 0. : 1.);
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

    // todo: freeze/unfreeze basis?

    /**
     * @brief branches by fixing a column
     */
    void branch() {
        info.n_iter_3cycle_current_node = 0;

        // gather branch node info
        const double obj_val_old = lp_solver_ptr->get_objective_value();
        const std::size_t j = reli_branch_ptr->get_branching_column();
        const double x_j_old = reli_branch_ptr->get_column_value(j);

        // branch
        ++info.n_branch_nodes;
        lp_solver_ptr->fix_column(j, reli_branch_ptr->get_up_cost(j) > reli_branch_ptr->get_down_cost(j));
        PACE_DEBUG_PRINTF("(branch)\n");
        stack.emplace(branch_node{j, obj_val_old, x_j_old, true});
    }

    /**
     * @brief given the optimal value of the current lp,
     * decides if we generate cutting planes,
     * or if we branch,
     * or if we backtrack.
     *
     * @return true optimal solution found
     * @return false otherwise
     */
    inline bool branch_and_bound_and_cut(std::vector<vertex_t> &ordering) {
        // test if the lp is optimal
        if (!lp_solver_ptr->is_optimal()) {
            return backtrack();
        }

        // test if we are worse than our current solution
        if (lp_solver_ptr->get_rounded_objective_value() >= upper_bound) {
            return backtrack();
        }

        // try to generate cutting planes
        bool cut_generated = true;
        if (lp_solver_ptr->is_integral() || info.n_iter_3cycle_current_node <= params.max_iter_3cycle_until_branch) {
            cut_generated = lp_solver_ptr->cut();
            params.delete_rows &= cut_generated;
            if (cut_generated) {
                if (params.delete_rows) lp_solver_ptr->delete_positive_slack_rows();
                return false;
            }
        }

        // at this point, we have an optimal solution to the ilp relaxation
        update_costs();

        // test if solution is integral, then we found a better solution
        if (lp_solver_ptr->is_integral() && !cut_generated) {
            build_ordering(ordering);
            lp_solver_ptr->fix_columns();
            return backtrack();
        }

        const crossing_number_t ub_old = upper_bound;
        if (heuristic.informed(*lp_solver_ptr, another_ordering) < ub_old) std::swap(ordering, another_ordering);
        info.relax_h_confidence = heuristic.get_confidence();

        // branch by fixing a column
        branch();
        return false;
    }

   public:
    /**
     * @brief solves the given instance exactly with a branch and cut algorithm
     */
    uint32_t operator()(std::vector<vertex_t> &ordering) {
        PACE_DEBUG_PRINTF("\nstart branch and cut\n");
        info.n_crossings_h = heuristic.uninformed(ordering);
        if (lower_bound() >= upper_bound) return upper_bound;

        if (!lp_solver_ptr) lp_solver_ptr = std::make_unique<highs_lp>(instance_);
        if (unsettled_pairs().size() == 0) {
            build_ordering(ordering);
            return upper_bound;
        }
        if (!reli_branch_ptr) reli_branch_ptr = std::make_unique<reliability_branching>(*lp_solver_ptr);

        const highs_lp_info &info_lp = lp_solver_ptr->get_info();
        PACE_DEBUG_PRINTF("start branch and cut\n");
        PACE_DEBUG_PRINTF("%11s=%11u, %11s=%11u, %11s=%11u\n",  //
                          "n cols", info_lp.n_cols,             //
                          "n initrows", info_lp.n_rows,         //
                          "n cand", info_lp.n_init_rows_candidates);

        // driver loop for branch and cut
        for (;;) {
            lp_solver_ptr->run();
            info.n_iter_3cycle_current_node += info_lp.new_3cycle_iter;
            if (branch_and_bound_and_cut(ordering)) break;
            info.depth = stack.size();
            PACE_DEBUG_PRINTF_INFO(info, info_lp);
            ++info.n_iterations;
        }

        PACE_DEBUG_PRINTF("end   branch and cut\n");
        return upper_bound;
    }

    //
    // getter
    //

    const branch_and_cut_info &get_info() {
        if (lp_solver_ptr) {
            info.n_rows = lp_solver_ptr->get_n_rows();
        }
        return info;
    }
};

};  // namespace pace

#endif