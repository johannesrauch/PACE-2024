#ifndef PACE_BRANCH_AND_CUT_HPP
#define PACE_BRANCH_AND_CUT_HPP

#include <glpk.h>
#include <math.h>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <functional>
#include <limits>
#include <stack>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "barycenter_heuristic.hpp"
#include "crossings.hpp"
#include "debug_printf.hpp"
#include "heuristics.hpp"
#include "instance.hpp"
#include "lp_wrapper.hpp"
#include "oracle.hpp"
#include "output.hpp"
#include "reliability_branching.hpp"
#include "topological_sort.hpp"

namespace pace {

struct branch_node {
    const std::size_t j;
    const double obj_val_old;
    const double x_j_old;
    bool fix_opposite;
};

struct branch_and_cut_info {
    const uint32_t &lower_bound;
    const uint32_t &upper_bound;

    std::size_t nof_iterations{0};
    std::size_t nof_branch_nodes{0};
    std::size_t nof_rows{0};

    uint32_t nof_crossings_h{0};
    uint32_t nof_crossings{0};
};

template <typename T, typename R>
class branch_and_cut {
   private:
    /// @brief the instance
    const pace::instance<T, R> &instance;

    /// @brief  n_free = graph.get_n_free()
    const std::size_t n_free;

    /// @brief lower bound on the optimal value of the instance
    const uint32_t &lower_bound;

    /// @brief current best upper bound on the optimal value of the instance
    uint32_t upper_bound{std::numeric_limits<uint32_t>::max()};

    double obj_val_old{0};

    /// @brief ordering corresponding to upper_bound
    std::vector<T> ordering;

    /**
     * @brief the restrictions graph for the instance.
     * the vertices of restriction_graph are the vertices of the free layer of the (instance) graph.
     * an arc (i, j) in restriction_graph means that i < j in the ordering.
     */
    digraph<T> restriction_graph;

    /**
     * @brief if magic[flat_index(u, v)] < 0, ~magic[flat_index(u, v)] is the value of x_uv in
     * an optimal solution. otherwise magic[flat_index(u, v)] is the variable index of x_uv.
     */
    std::vector<int> magic;

    /// @brief contains pairs not settled by the oracle
    std::vector<std::pair<T, T>> unsettled;

    /// @brief for a dfs search tree
    std::stack<branch_node> stack;

    /// @brief lp wrapper
    highs_wrapper<T> *lp_solver{nullptr};

    reliability_branching<T> *reli_branch{nullptr};

    branch_and_cut_info info;

   public:
    /**
     * @brief constructs and initializes the branch and cut solver
     */
    branch_and_cut(const pace::instance<T, R> &instance)
        : instance(instance),
          n_free(instance.graph().get_n_free()),
          lower_bound(instance.get_lower_bound()),
          ordering(n_free),
          restriction_graph(n_free),
          info{lower_bound, upper_bound} {
        assert(n_free > 0);
    }

    ~branch_and_cut() {
        if (lp_solver != nullptr) {
            delete lp_solver;
            lp_solver = nullptr;
        }
        if (reli_branch != nullptr) {
            delete reli_branch;
            reli_branch = nullptr;
        }
    }

    // delete copy constructor, move constructor, copy assignment and move assignment
    branch_and_cut(const branch_and_cut<T, R> &rhs) = delete;
    branch_and_cut(branch_and_cut<T, R> &&rhs) = delete;
    branch_and_cut<T, R> &operator=(const branch_and_cut<T, R> &rhs) = delete;
    branch_and_cut<T, R> &operator=(branch_and_cut<T, R> &&rhs) = delete;

   private:
    /**
     * @brief constructs the restrictions graph into restriction_graph.
     * the topological sort of `restriction_graph` gives an ordering of the vertices of the free
     * layer of `graph`, which we store in `ordering`.
     *
     * @pre lp solution integral
     */
    void compute_ordering() {
        assert(lp_solver->is_integral());
        restriction_graph.rollback();

        // build the constraint graph
        for (std::size_t i = 0; i < unsettled.size(); ++i) {
            const double x = lp_solver->get_column_value(i);
            const auto &[u, v] = unsettled[i];
            if (x < 0.5) {
                restriction_graph.add_arc(v, u);
            } else {
                restriction_graph.add_arc(u, v);
            }
        }

        // the ordering computed by the lp is the topological sort
        bool acyclic = topological_sort(restriction_graph, ordering);
        (void)acyclic;  // suppress unused warning
        assert(acyclic);
        uint32_t new_upper_bound = lp_solver->get_rounded_objective_value();
        assert(new_upper_bound < upper_bound);
        assert(new_upper_bound == number_of_crossings(instance.graph(), ordering));

        // try to improve new solution
        new_upper_bound = shift_heuristic{instance}(ordering, new_upper_bound);
        upper_bound = new_upper_bound;
    }

    /// @brief updates pseudo-costs
    void update_costs() {
        if (stack.empty()) return;
        const auto &[j, obj_val_old, x_j_old, _] = stack.top();
        const double obj_val_new = lp_solver->get_objective_value();
        if (lp_solver->get_column_value(j) > 0.5) {
            reli_branch->update_up_cost(j, obj_val_new, obj_val_old, x_j_old);
        } else {
            reli_branch->update_down_cost(j, obj_val_new, obj_val_old, x_j_old);
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
                lp_solver->fix_column(node.j, lp_solver->get_column_value(node.j) > 0.5 ? 0. : 1.);
                PACE_DEBUG_PRINTF("(backtrack)\n");
                node.fix_opposite = false;
                return false;
            } else {
                lp_solver->unfix_column(node.j);
                PACE_DEBUG_PRINTF("unfixed variable %5d\n(backtrack)\n", node.j);
                stack.pop();
            }
        }
        return true;
    }

    /// @brief branches by fixing a column
    void branch() {
        // gather branch node info
        const double obj_val_old = lp_solver->get_objective_value();
        const std::size_t j = reli_branch->get_branching_column();
        const double x_j_old = reli_branch->get_column_value(j);

        // branch
        ++info.nof_branch_nodes;
        lp_solver->fix_column(j, reli_branch->get_up_cost(j) > reli_branch->get_down_cost(j));
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
        if (!lp_solver->is_optimal()) {
            return backtrack();
        }

        // test if we are worse than our current solution
        const double value = lp_solver->get_objective_value();
        if (static_cast<double>(upper_bound) <= value) {
            return backtrack();
        }

        // try to generate cutting planes
        bool successful = lp_solver->cut();
        if (successful) {
            if (stack.empty()) lp_solver->delete_positive_slack_rows();
            return false;
        }

        // at this point, we have an optimal solution to the ilp relaxation
        update_costs();

        // test if solution is integral, then we found a better solution
        if (lp_solver->is_integral()) {
            compute_ordering();
            lp_solver->fix_columns(upper_bound);
            return backtrack();
        }

        // branch by fixing a column
        branch();
        return false;
    }

   public:
    /**
     * @brief solves the given instance exactly with a branch and cut algorithm
     *
     * @param do_print to print or not to print (the solution)
     */
    template <bool DO_OUTPUT_PRINT = true>
    uint32_t operator()() {
        PACE_DEBUG_PRINTF("\n\nstart heuristic\n");
        upper_bound = heuristics(instance, ordering);
        info.nof_crossings_h = upper_bound;
        PACE_DEBUG_PRINTF("end   heuristic\n");
        if (lower_bound >= upper_bound) return upper_bound;

        if (lp_solver == nullptr) {
            PACE_DEBUG_PRINTF("start oracle\n");
            oracle<T, R> oracle(instance, upper_bound);
            const uint32_t objective_offset = oracle.build(restriction_graph, magic, unsettled);
            PACE_DEBUG_PRINTF("end   oracle\n");

            lp_solver = new highs_wrapper<T>(instance, magic, unsettled, objective_offset, upper_bound);
            reli_branch = new reliability_branching<T>(*lp_solver);
        }

        // driver loop for branch and cut
        const highs_wrapper_info<T> &info_lp = lp_solver->get_info();
        (void) info_lp;
        PACE_DEBUG_PRINTF("start branch and cut\n");
        PACE_DEBUG_PRINTF("%11s=%11u\n", "nof cols", unsettled.size());
        do {
            lp_solver->run();
            PACE_DEBUG_INFO(info, info_lp);
            ++info.nof_iterations;
        } while (!branch_and_bound_and_cut());
        PACE_DEBUG_PRINTF("end   branch and cut\n");

        PACE_DEBUG_PRINTF("OPT: %u\n", upper_bound);
        if constexpr (DO_OUTPUT_PRINT) {
            print_output(instance.graph().get_n_fixed(), ordering);
        }
        return upper_bound;
    }

    //
    // getter
    //

    /**
     * @brief get number of crossings of current ordering
     */
    uint32_t get_nof_crossings() { return upper_bound; }

    /**
     * @brief get current best ordering
     *
     * @return const std::vector<T>&
     */
    const std::vector<T> &get_ordering() const { return ordering; }

    const branch_and_cut_info &get_info() {
        if (lp_solver != nullptr) {
            info.nof_rows = lp_solver->get_nof_rows();
        }
        info.nof_crossings = upper_bound;
        return info;
    }
};

};  // namespace pace

#endif