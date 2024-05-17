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
    std::size_t max_nodes{std::numeric_limits<std::size_t>::max()};
    std::size_t n_nodes_until_rins{25};
    bool do_rins{true};
    bool do_uninformed_h{true};
    bool do_uninformed_rins{false};
    bool do_lsearch{false};
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

    heuristics heuristic;

   public:
    /**
     * @brief constructs and initializes the branch and cut solver
     */
    branch_and_cut(instance &instance_,
                   branch_and_cut_params params = branch_and_cut_params())
        : instance_view(instance_),
          info{lower_bound(), upper_bound},
          params(params),
          heuristic(instance_) {
        assert(n_free > 0);
    }

    // delete copy constructor, move constructor, copy assignment and move
    // assignment
    branch_and_cut(const branch_and_cut &rhs) = delete;
    branch_and_cut(branch_and_cut &&rhs) = delete;
    branch_and_cut &operator=(const branch_and_cut &rhs) = delete;
    branch_and_cut &operator=(branch_and_cut &&rhs) = delete;

   private:
    /**
     * @brief constructs the restrictions graph into restriction_graph.
     * the topological sort of `restriction_graph` gives an ordering of the
     * vertices of the free layer of `graph`, which we store in `ordering`.
     */
    void build_ordering(std::vector<vertex_t> &ordering);

    /**
     * @brief updates pseudo-costs
     */
    void update_costs();

    //
    // branch and bound and cut methods
    //

    /**
     * @brief if there's a variable on the stack,
     * the function pops it, and we fix the opposite value or we backtrack
     * again.
     *
     * @return true if stack is empty
     * @return false otherwise
     */
    bool backtrack();

    /**
     * @brief branches by fixing a column
     */
    bool branch();

    /**
     * @brief given the optimal value of the current lp,
     * decides if we generate cutting planes,
     * or if we branch,
     * or if we backtrack.
     *
     * @return true optimal solution found
     * @return false otherwise
     */
    bool branch_and_bound_and_cut(std::vector<vertex_t> &ordering);

   public:
    /**
     * @brief solves the given instance exactly with a branch and cut algorithm
     */
    crossing_number_t operator()(std::vector<vertex_t> &ordering);

    //
    // getter
    //

    const branch_and_cut_info &get_info() { return info; }
};

};  // namespace pace

#endif