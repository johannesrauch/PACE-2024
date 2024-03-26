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
#include "topological_sort.hpp"

namespace pace {

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

    /// @brief for branching
    std::stack<std::pair<int, bool>> stack;

    /// @brief lp wrapper
    highs_wrapper<T> *lp_solver{nullptr};

   public:
    /**
     * @brief constructs and initializes the branch and cut solver
     */
    branch_and_cut(const pace::instance<T, R> &instance)
        : instance(instance),
          n_free(instance.graph().get_n_free()),
          lower_bound(instance.get_lower_bound()),
          ordering(n_free),
          restriction_graph(n_free) {
        assert(n_free > 0);
    }

    ~branch_and_cut() {
        if (lp_solver != nullptr) {
            delete lp_solver;
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
        assert(lp_solver->is_integral() == -1);
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
        const uint32_t new_upper_bound = lp_solver->get_rounded_objective_value();
        assert(new_upper_bound < upper_bound);
        assert(new_upper_bound == number_of_crossings(instance.graph(), ordering));
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
                lp_solver->fix_column(j, lp_solver->get_variable_value(j) > 0.5 ? 0. : 1.);
                PACE_DEBUG_PRINTF("(backtrack)\n");
                stack.emplace(j, false);
                break;
            } else {
                lp_solver->unfix_column(j);
                PACE_DEBUG_PRINTF("unfixed variable %5d\n(backtrack)\n", j);
            }
        }
        return !fix_opposite;
    }

    /// @brief branches by fixing a column
    void branch(const int &j) {
        // todo: more sophisticated fixing
        lp_solver->fix_column(j, 0.);
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

        const int j = lp_solver->is_integral();
        // test if solution is integral, then we found a better solution
        if (j == -1) {
            compute_ordering();
            PACE_DEBUG_PRINTF("\tfound a better solution with %u crossings\n", upper_bound);
            lp_solver->fix_columns(upper_bound);
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
    template <bool DO_OUTPUT_PRINT = true>
    uint32_t operator()() {
        PACE_DEBUG_PRINTF("start heuristic\n");
        upper_bound = heuristics(instance, ordering);
        PACE_DEBUG_PRINTF("end   heuristic\n");
        if (lower_bound >= upper_bound) return upper_bound;

        if (lp_solver == nullptr) {
            PACE_DEBUG_PRINTF("start oracle\n");
            oracle<T, R> oracle(instance, upper_bound);
            const uint32_t objective_offset = oracle.build(restriction_graph, magic, unsettled);
            PACE_DEBUG_PRINTF("end   oracle\n");

            lp_solver = new highs_wrapper<T>(instance, magic, objective_offset);
        }

        // driver loop for branch and cut
        std::size_t iteration = 0;
        do {
            (void)iteration;
            PACE_DEBUG_PRINTF(
                "iteration=%5u, "
                "depth=%5u, "
                "objective_val=%5.0f, "
                "upper_bound=%5u, "
                "nof_rows=%5u\n",
                iteration++,                       //
                stack.size(),                      //
                lp_solver->get_objective_value(),  //
                upper_bound,                       //
                lp_solver->get_nof_rows());
            lp_solver->run();
        } while (!branch_and_bound_and_cut());

        PACE_DEBUG_PRINTF("OPTIMAL VALUE:   %u\n", upper_bound);
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

    std::size_t get_nof_rows() {
        if (lp_solver != nullptr) {
            return lp_solver->get_nof_rows();
        } else {
            return 0;
        }
    }
};

};  // namespace pace

#endif