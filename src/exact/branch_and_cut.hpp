#ifndef PACE_EXACT_BRANCH_AND_CUT_HPP
#define PACE_EXACT_BRANCH_AND_CUT_HPP

#include <cfloat>
#include <cmath>
#include <functional>
#include <stack>
#include <unordered_map>

#include "exact/lp_wrapper.hpp"
#include "exact/reliability_branching.hpp"
#include "heuristic/heuristics.hpp"
#include "io/output.hpp"
#include "log/debug_printf.hpp"
#include "model/instance.hpp"
#include "utils/crossings_utils.hpp"
#include "utils/topological_sort.hpp"

namespace pace {

class branch_and_cut : public instance_view {
    /**
     * @brief ordering corresponding to upper_bound
     */
    std::vector<vertex_t> ordering;

    struct branch_node {
        const std::size_t j;
        const double obj_val_old;
        const double x_j_old;
        bool fix_opposite;
    };

    /**
     * @brief for a dfs through search tree
     */
    std::stack<branch_node> stack;

    /**
     * @brief interface to the ilp relaxation solver
     */
    std::unique_ptr<highs_wrapper> lp_solver_ptr;

    std::unique_ptr<reliability_branching> reli_branch_ptr;

    branch_and_cut_info info;

   public:
    /**
     * @brief constructs and initializes the branch and cut solver
     */
    branch_and_cut(instance &instance_) : instance_view(instance_), ordering(n_free), info{lower_bound(), upper_bound} {
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
     *
     * @pre lp solution integral
     */
    void compute_ordering() {
        assert(lp_solver_ptr->is_integral());
        digraph &restr_graph{restriction_graph()};
        restr_graph.rollback();

        // build the constraint graph
        std::size_t i = 0;
        for (const auto &[u, v] : unsettled_pairs()) {
            const double x = lp_solver_ptr->get_column_value(i);
            if (x < 0.5) {
                restr_graph.add_arc(v, u);
            } else {
                restr_graph.add_arc(u, v);
            }
            ++i;
        }

        // the ordering computed by the lp is the topological sort
        bool acyclic = topological_sort(restr_graph, ordering);
        (void)acyclic;  // suppress unused warning
        assert(acyclic);

        crossing_number_t n_crossings = lp_solver_ptr->get_rounded_objective_value();
        assert(n_crossings < upper_bound);
        assert(n_crossings == number_of_crossings(graph, ordering));

        // try to improve new solution
        n_crossings = shift_heuristic{instance_}(ordering, n_crossings);
        update_upper_bound(n_crossings);
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
        while (!stack.empty()) {
            branch_node &node = stack.top();
            if (node.fix_opposite) {
                lp_solver_ptr->fix_column(node.j, lp_solver_ptr->get_column_value(node.j) > 0.5 ? 0. : 1.);
                PACE_DEBUG_PRINTF("(backtrack)\n");
                node.fix_opposite = false;
                return false;
            } else {
                lp_solver_ptr->unfix_column(node.j);
                PACE_DEBUG_PRINTF("unfixed variable %5d\n(backtrack)\n", node.j);
                stack.pop();
            }
        }
        return true;
    }

    /**
     * @brief branches by fixing a column
     */
    void branch() {
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
    inline bool branch_and_bound_and_cut() {
        // test if the lp is optimal
        if (!lp_solver_ptr->is_optimal()) {
            return backtrack();
        }

        // test if we are worse than our current solution
        const double value = lp_solver_ptr->get_objective_value();
        if (static_cast<double>(upper_bound) <= value) {
            return backtrack();
        }

        // try to generate cutting planes
        bool successful = lp_solver_ptr->cut();
        if (successful) {
            if (stack.empty()) lp_solver_ptr->delete_positive_slack_rows();
            return false;
        }

        // at this point, we have an optimal solution to the ilp relaxation
        update_costs();

        // test if solution is integral, then we found a better solution
        if (lp_solver_ptr->is_integral()) {
            compute_ordering();
            lp_solver_ptr->fix_columns();
            return backtrack();
        }

        // branch by fixing a column
        branch();
        return false;
    }

   public:
    /**
     * @brief solves the given instance exactly with a branch and cut algorithm
     */
    uint32_t operator()() {
        info.n_crossings_h = heuristics(instance_, ordering);
        if (lower_bound() >= upper_bound) return upper_bound;

        if (!lp_solver_ptr) lp_solver_ptr = std::make_unique<highs_wrapper>(instance_);
        if (!reli_branch_ptr) reli_branch_ptr = std::make_unique<reliability_branching>(*lp_solver_ptr);

#ifndef NDEBUG
        const highs_wrapper_info &info_lp = lp_solver_ptr->get_info();
        (void)info_lp;
        PACE_DEBUG_PRINTF("start branch and cut\n");
        PACE_DEBUG_PRINTF("%11s=%11u, %11s=%11u\n",               //
                          "nof cols", lp_solver->get_nof_cols(),  //
                          "nof i rows", lp_solver->get_nof_rows());
#endif

        // driver loop for branch and cut
        do {
            lp_solver_ptr->run();
            PACE_DEBUG_INFO(info, info_lp);
            ++info.n_iterations;
        } while (!branch_and_bound_and_cut());
        PACE_DEBUG_PRINTF("end   branch and cut\n");

        return upper_bound;
    }

    //
    // getter
    //

    const std::vector<vertex_t> &get_ordering() const { return ordering; }

    const branch_and_cut_info &get_info() {
        if (lp_solver_ptr) {
            info.n_rows = lp_solver_ptr->get_nof_rows();
        }
        return info;
    }
};

};  // namespace pace

#endif