#ifndef PACE_BRANCH_AND_CUT_HPP
#define PACE_BRANCH_AND_CUT_HPP

#include <glpk.h>
#include <math.h>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <functional>
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
#include "topological_sort.hpp"

namespace pace {

template <typename T, typename R>
class branch_and_cut {
   private:
    /// @brief a bipartite graph that resembles an instance of one-sided crossing minimization
    const bipartite_graph<T> &graph;

    /// @brief  n_free = graph.get_n_free()
    const std::size_t n_free;

    /// @brief lp wrapper
    highs_wrapper<T> lp_solver;

    /// @brief lower bound on the optimal value of the instance
    const uint32_t &lower_bound;

    /// @brief current best upper bound on the optimal value of the instance
    uint32_t upper_bound{0};

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

    std::vector<std::pair<T, T>> unsettled;

    /// @brief for branching
    std::stack<std::pair<int, bool>> stack;

   public:
    /**
     * @brief constructs and initializes the branch and cut solver     *
     *
     * @param graph input graph
     */
    branch_and_cut(const instance<T, R> &instance)
        : graph(instance.graph()),
          n_free(graph.get_n_free()),
          lp_solver(instance, magic),
          lower_bound(instance.get_lower_bound()),
          ordering(n_free),
          restriction_graph(n_free) {
        assert(n_free > 0);

        PACE_DEBUG_PRINTF("start heuristic\n");
        upper_bound = heuristics(instance, ordering);
        PACE_DEBUG_PRINTF("end   heuristic\n");

        if (lower_bound >= upper_bound) {
            PACE_DEBUG_PRINTF("start oracle\n");
            oracle<T, R> oracle(instance, upper_bound);
            oracle.build(restriction_graph, magic, unsettled);
            PACE_DEBUG_PRINTF("end   oracle\n");
        }
    }

    // delete copy constructor, move constructor, copy assignment and move assignment
    branch_and_cut(const branch_and_cut<T, R> &rhs) = delete;
    branch_and_cut(branch_and_cut<T, R> &&rhs) = delete;
    branch_and_cut<T, R> &operator=(const branch_and_cut<T, R> &rhs) = delete;
    branch_and_cut<T, R> &operator=(branch_and_cut<T, R> &&rhs) = delete;

   private:
    //
    // solution related methods
    //

    /**
     * @brief constructs the restrictions graph into restriction_graph.
     * the topological sort of `restriction_graph` gives an ordering of the vertices of the free
     * layer of `graph`, which we store in `ordering`.
     *
     * @pre lp solution integral
     */
    void compute_ordering() {
        assert(lp_solver.is_integral() == -1);
        restriction_graph.rollback();

        // build the constraint graph
        for (std::size_t i = 0; i < unsettled.size(); ++i) {
            const double x = lp_solver.get_column_value(i);
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
        const uint32_t new_upper_bound = lp_solver.get_rounded_objective_value();
        assert(new_upper_bound < upper_bound);
        upper_bound = new_upper_bound;
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
        int j{0};
        bool fix_opposite = false;
        for (; !stack.empty();) {
            std::tie(j, fix_opposite) = stack.top();
            stack.pop();

            if (fix_opposite) {
                lp_solver.fix_column(j, lp_solver.get_variable_value(j) > 0.5 ? 0. : 1.);
                PACE_DEBUG_PRINTF("(backtrack)\n");
                stack.emplace(j, false);
                break;
            } else {
                lp_solver.unfix_column(j);
                PACE_DEBUG_PRINTF("unfixed variable %5d\n(backtrack)\n", j);
            }
        }
        return !fix_opposite;
    }

    /// @brief branches by fixing a column
    void branch(const int &j) {
        // todo: more sophisticated fixing
        lp_solver.fix_column(j, 0.);
        PACE_DEBUG_PRINTF("(branch)\n");
        stack.emplace(j, true);
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
        if (!lp_solver.is_optimal()) {
            return backtrack();
        }

        // test if we are worse than our current solution
        const double value = lp_solver.get_objective_value();
        if (static_cast<double>(upper_bound) <= value) {
            return backtrack();
        }

        // try to generate cutting planes
        bool successful = lp_solver.cut();
        if (successful) {
            if (stack.empty()) lp_solver.delete_positive_slack_rows();
            return false;
        }

        const int j = lp_solver.is_integral();
        // test if solution is integral, then we found a better solution
        if (j == -1) {
            compute_ordering();
            lp_solver.fix_columns(upper_bound);
            return backtrack();
        }

        // todo: improve solution with heuristic
        // branch by fixing a column
        branch(j);
        return false;
    }

   public:
    /**
     * @brief solves the given instance exactly with a branch and cut algorithm
     *
     * @param do_print to print or not to print (the solution)
     */
    void run(bool do_print = true) {
        if (upper_bound == 0) return;

        // solve lp
        lp_solver.run();
        assert(lp_solver.is_optimal());

        // get value of current optimal solution and call branch_and_bound_and_cut with it
        bool is_optimum = branch_and_bound_and_cut();

        // driver loop
        for (std::size_t iteration = 2; !is_optimum; ++iteration) {
            PACE_DEBUG_PRINTF("iteration=%5llu, depth=%5llu, upper_bound=%5lu, nof_rows=%5d\n",
                              iteration, stack.size(), upper_bound, lp_solver.get_nof_rows());

            // solve lp
            lp_solver.run();

            // branch or cut
            is_optimum = branch_and_bound_and_cut();
        }

        // output
        PACE_DEBUG_PRINTF("OPTIMAL VALUE: %lu\n", upper_bound);
        if (do_print) {
            print_output(graph.get_n_fixed(), ordering);
        }
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
};

};  // namespace pace

#endif