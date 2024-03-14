#ifndef PACE2024_BRANCH_AND_CUT_HPP
#define PACE2024_BRANCH_AND_CUT_HPP

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
#include "bipartite_graph.hpp"
#include "crossings.hpp"
#include "debug_printf.hpp"
#include "lp_wrapper.hpp"
#include "matrix.hpp"
#include "median_heuristic.hpp"
#include "output.hpp"
#include "topological_sort.hpp"

namespace pace2024 {

template <typename T>
class branch_and_cut {
   private:
    /// @brief a bipartite graph that resembles an instance of one-sided crossing minimization
    const bipartite_graph<T> &graph;

    /// @brief  n1 = graph.get_n1()
    const std::size_t n1;

    /// @brief pointer to lp solver
    highs_wrapper lp_solver;

    /// @brief best upper bound on the optimal value of the instance
    uint32_t upper_bound{0};

    /// @brief ordering corresponding to upper_bound
    std::vector<T> ordering;

    /**
     * @brief the restrictions graph for the instance.
     * the vertices of restriction_graph are the vertices of the free layer of the (instance) graph.
     * an arc (i, j) in restriction_graph means that i < j in the ordering.
     */
    digraph<T> restriction_graph;

    /// @brief for branching
    std::stack<std::pair<int, bool>> stack;

   public:
    /**
     * @brief constructs and initializes the branch and cut solver     *
     *
     * @param graph input graph
     */
    branch_and_cut(bipartite_graph<T> &graph)
        : graph(graph),
          n1(graph.get_n_free()),
          lp_solver(graph),
          ordering(n1),
          restriction_graph(n1) {
        assert(n1 > 0);
    }

    // delete copy constructor, move constructor, copy assignment and move assignment
    branch_and_cut(const branch_and_cut &rhs) = delete;
    branch_and_cut(branch_and_cut &&rhs) = delete;
    branch_and_cut &operator=(const branch_and_cut &rhs) = delete;
    branch_and_cut &operator=(branch_and_cut &&rhs) = delete;

   private:
    //
    // solution related methods
    //

    /**
     * @brief constructs the restrictions graph into restriction_graph.
     * the topological sort of `restriction_graph` gives an ordering of the vertices of the free layer of `graph`,
     * which we store in `ordering`.
     *
     * @pre lp solution integral
     */
    void compute_ordering() {
        assert(lp_solver.is_integral() == -1);
        restriction_graph.clear_arcs();

        // build the constraint graph
        int k = 0;
        for (T i = 0; i < n1; ++i) {
            for (T j = i + 1; j < n1; ++j) {
                const double x = lp_solver.get_variable_value(k);
                if (x < 0.5) {
                    restriction_graph.add_arc(j, i);
                } else {
                    restriction_graph.add_arc(i, j);
                }
                ++k;
            }
        }

        // the ordering computed by the lp is the topological sort
        bool acyclic = topological_sort(restriction_graph, ordering);
        (void)acyclic;  // suppress unused warning
        assert(acyclic);
        const uint32_t new_upper_bound = static_cast<uint32_t>(lp_solver.get_rounded_objective_value());
        assert(new_upper_bound < upper_bound);
        upper_bound = new_upper_bound;
    }

    //
    // heuristic related methods
    //

    /**
     * @brief runs heuristics to get an initial value for `upper_bound` and `ordering`
     */
    void run_heuristics() {
        PACE2024_DEBUG_PRINTF("start heuristic\n");

        median_heuristic(graph, ordering).run();
        upper_bound = number_of_crossings(graph, ordering);

        std::vector<T> ordering_;
        barycenter_heuristic(graph, ordering_).run();
        const uint32_t upper_bound_ = number_of_crossings(graph, ordering_);

        if (upper_bound_ < upper_bound) {
            upper_bound = upper_bound_;
            ordering = ordering_;
        }

        PACE2024_DEBUG_PRINTF("end   heuristic\n");
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
                PACE2024_DEBUG_PRINTF("(backtrack)\n");
                stack.emplace(j, false);
                break;
            } else {
                lp_solver.unfix_column(j);
                PACE2024_DEBUG_PRINTF("unfixed variable %5d\n(backtrack)\n", j);
            }
        }
        return !fix_opposite;
    }

    /// @brief branches by fixing a column
    void branch(const int &j) {
        // todo: more sophisticated fixing
        lp_solver.fix_column(j, 0.);
        PACE2024_DEBUG_PRINTF("(branch)\n");
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
            lp_solver.fix_columns(static_cast<int>(upper_bound));
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
    void solve(bool do_print = true) {
        run_heuristics();  // to initialize `upper_bound` and `ordering`
        if (upper_bound == 0) return;

        // solve lp
        lp_solver.solve();
        assert(lp_solver.is_optimal());

        // get value of current optimal solution and call branch_and_bound_and_cut with it
        bool is_optimum = branch_and_bound_and_cut();

        // driver loop
        for (std::size_t iteration = 2; !is_optimum; ++iteration) {
            PACE2024_DEBUG_PRINTF("iteration=%5llu, depth=%5llu, upper_bound=%5lu, nof_rows=%5d\n",
                                  iteration,
                                  stack.size(),
                                  upper_bound,
                                  lp_solver.get_nof_rows());

            // solve lp
            lp_solver.solve();

            // branch or cut
            is_optimum = branch_and_bound_and_cut();
        }

        // output
        PACE2024_DEBUG_PRINTF("OPTIMAL VALUE: %lu\n", upper_bound);
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
    uint32_t get_nof_crossings() {
        return upper_bound;
    }

    /**
     * @brief get current best ordering
     *
     * @return const std::vector<T>
     */
    const std::vector<T> get_ordering() const {
        return ordering;
    }
};

};  // namespace pace2024

#endif