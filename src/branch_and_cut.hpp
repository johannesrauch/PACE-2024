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

#include "bipartite_graph.hpp"
#include "crossing_number.hpp"
#include "debug_printf.hpp"
#include "lp_wrapper.hpp"
#include "matrix.hpp"
#include "median_heuristic.hpp"
#include "output.hpp"
#include "topological_sort.hpp"

namespace pace2024 {

template <typename T, typename R>
class branch_and_cut {
   private:
    /// @brief a bipartite graph that resembles an instance of one-sided crossing minimization
    const general_bipartite_graph<T> &graph;

    /// @brief  n1 = graph.get_n1()
    const std::size_t n1;

    /// @brief pointer to lp solver
    lp_wrapper *const lp_solver;

    /// @brief best upper bound on the optimal value of the instance
    R upper_bound{0};

    /// @brief ordering corresponding to upper_bound
    std::vector<T> ordering;

    /**
     * @brief the restrictions graph for the instance.
     * the vertices of digraph are the vertices of the free layer of the (instance) graph.
     * an arc (i, j) in digraph means that i < j in the ordering.
     */
    general_digraph<T> digraph;

    /// @brief for branching
    std::stack<std::pair<int, bool>> stack;

   public:
    /**
     * @brief constructs and initializes the branch and cut solver     *
     *
     * @param graph input graph
     */
    branch_and_cut(general_bipartite_graph<T> &graph)
        : graph(graph),
          n1(graph.get_n1()),
          lp_solver(new clp_wrapper(graph)),
          ordering(n1),
          digraph(n1) {
        assert(n1 > 0);
        run_heuristics();  // to initialize `upper_bound` and `ordering`
    }

    // delete copy constructor, move constructor, copy assignment and move assignment
    branch_and_cut(const branch_and_cut &rhs) = delete;
    branch_and_cut(branch_and_cut &&rhs) = delete;
    branch_and_cut &operator=(const branch_and_cut &rhs) = delete;
    branch_and_cut &operator=(branch_and_cut &&rhs) = delete;

    ~branch_and_cut() {
        delete lp_solver;
    }

   private:
    //
    // solution related methods
    //

    /**
     * @brief constructs the restrictions graph into digraph.
     * the topological sort of `digraph` gives an ordering of the vertices of the free layer of `graph`,
     * which we store in `ordering`.
     *
     * @pre lp solution integral
     */
    void compute_ordering() {
        assert(lp_solver->is_integral() == 0);
        digraph.delete_arcs();

        // build the constraint graph
        int k = 1;
        for (T i = 0; i < n1; ++i) {
            for (T j = i + 1; j < n1; ++j) {
                const double x = lp_solver->get_column_value(k);
                if (x < 0.5) {
                    digraph.add_arc(j, i);
                } else {
                    digraph.add_arc(i, j);
                }
                ++k;
            }
        }

        // the ordering computed by the lp is the topological sort
        bool acyclic = topological_sort(digraph, ordering);
        (void)acyclic;  // suppress unused warning
        assert(acyclic);
        upper_bound = static_cast<R>(lp_solver->get_rounded_objective_value());
    }

    //
    // heuristic related methods
    //

    /**
     * @brief runs heuristics to get an initial value for `upper_bound` and `ordering`
     */
    void run_heuristics() {
        PACE2024_DEBUG_PRINTF("start heuristic\n");
        upper_bound = probabilistic_median_heuristic<T, R>(graph, ordering).run();
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
                lp_solver->fix_column(j, lp_solver->get_column_value(j) > 0.5 ? 0. : 1.);
                PACE2024_DEBUG_PRINTF("(backtrack)\n");
                stack.emplace(j, false);
                break;
            } else {
                lp_solver->unfix_column(j);
                PACE2024_DEBUG_PRINTF("unfixed variable %5d\n(backtrack)\n", j);
            }
        }
        return !fix_opposite;
    }

    /// @brief branches by fixing a column
    void branch(const int &j) {
        // todo: more sophisticated fixing
        lp_solver->fix_column(j, 0.);
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
        if (successful) return false;

        const int j = lp_solver->is_integral();
        // test if solution is integral, then we found a better solution
        if (j == 0) {
            compute_ordering();
            lp_solver->fix_columns(static_cast<int>(upper_bound));
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
        // solve lp
        lp_solver->solve(false);
        assert(lp_solver->is_optimal());

        // get value of current optimal solution and call branch_and_bound_and_cut with it
        bool is_optimum = branch_and_bound_and_cut();

        // driver loop
        for (std::size_t iteration = 0; !is_optimum; ++iteration) {
            PACE2024_DEBUG_PRINTF("depth=%llu, upper_bound=%lu, nof_rows=%d\n",
                                  stack.size(),
                                  upper_bound,
                                  lp_solver->get_nof_rows());

            // solve lp
            bool delete_rows_after = stack.size() == 0;
            lp_solver->solve(delete_rows_after);

            // branch or cut
            is_optimum = branch_and_bound_and_cut();
        }

        // output
        PACE2024_DEBUG_PRINTF("OPTIMAL VALUE: %lu\n", upper_bound);
        if (do_print) {
            print_output(graph.get_n0(), ordering);
        }
    }

    //
    // getter
    //

    /**
     * @brief get number of crossings of current ordering
     *
     * @return R
     */
    R get_nof_crossings() {
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