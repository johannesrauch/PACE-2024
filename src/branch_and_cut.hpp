#ifndef PACE2024_BRANCH_AND_CUT_HPP
#define PACE2024_BRANCH_AND_CUT_HPP

#ifndef PACE2024_CONST_NOF_CYCLE_CONSTRAINTS
#define PACE2024_CONST_NOF_CYCLE_CONSTRAINTS 128
#endif

// 1e-7 is the default tolerance in glpk for feasibility
#ifndef PACE2024_CONST_EPSILON
#define PACE2024_CONST_EPSILON 1e-7
#endif

#ifndef PACE2024_CONST_NOF_BUCKETS
#define PACE2024_CONST_NOF_BUCKETS 1000
#endif

#include <glpk.h>
#include <math.h>

#include <functional>
#include <stack>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "bipartite_graph.hpp"
#include "crossing_number.hpp"
#include "debug_printf.hpp"
#include "matrix.hpp"
#include "median_heuristic.hpp"
#include "output.hpp"
#include "topological_sort.hpp"

namespace pace2024 {

template <typename T, typename R>
class branch_and_cut {
   private:
    //
    // instance related attributes
    //

    /**
     * @brief a bipartite graph that resembles an instance of
     * one-sided crossing minimization
     */
    const general_bipartite_graph<T> &graph;

    /**
     * @brief number of vertices in free layer of the graph
     */
    const std::size_t n1;

    /**
     * @brief n1_choose_2 = n1*(n1-1)/2
     */
    const std::size_t n1_choose_2;

    //
    // lp related attributes
    //

    /**
     * @brief const pointer to (mutable) lp solver
     */
    glp_prob *const lp;

    /**
     * @brief parameters for glp_simplex
     */
    glp_smcp params;

    /**
     * @brief obj_val_offset + the obj_val of the lp = number of crossings
     */
    R obj_val_offset{0};

    /**
     * @brief stores an optimal solution to the lp without constraints
     * used for permanent variable fixing
     */
    std::vector<uint8_t> initial_solution;

    //
    // solution related attributes
    //

    /**
     * @brief lower bound for optimal value of the instance
     * given by computing sum min(c_ij, c_ji)
     */
    R lower_bound{0};

    /**
     * @brief upper bound for optimal value of the instance, that is,
     * number of crossings of current best solution.
     * - initialized through heuristic
     * - updated over time through this solver
     */
    R upper_bound{0};

    /**
     * @brief ordering achieving the current best number of crossings
     */
    std::vector<T> ordering;

    /**
     * @brief the restrictions graph for the instance.
     * the vertices of digraph are the vertices of the free layer of the
     * (instance) graph.
     * an arc (i, j) in digraph means that i < j in the ordering.
     */
    general_digraph<T> digraph;

    //
    // solver related attributes
    //

    /**
     * @brief stack used for branch and cut
     */
    std::stack<int> stack;

    /**
     * @brief buckets for bucket sorting violated 3-cycle ieqs
     */
    std::vector<std::tuple<T, T, T, bool>> buckets{PACE2024_CONST_NOF_BUCKETS};

   public:
    /**
     * @brief constructs and initializes the branch and cut solver:
     * - adds columns (variables)
     * - adds hypercube bounds
     * - fixes some variables, too
     * - computes a preliminary solution with median heuristic
     *
     * @param instance input
     */
    branch_and_cut(general_bipartite_graph<T> &graph, int msg_level = GLP_MSG_OFF)
        : graph(graph),
          n1(graph.get_n1()),
          n1_choose_2(n1 * (n1 - 1) / 2),
          lp(glp_create_prob()),
          initial_solution(n1_choose_2 + 1),
          ordering(n1),
          digraph(n1) {
        assert(n1 > 0);

        // initialize glp_simplex parameters and set lp objective to minimization
        glp_init_smcp(&params);
        params.msg_lev = msg_level;
        glp_set_obj_dir(lp, GLP_MIN);

        // compute the crossing numbers and add respective variables to the lp
        PACE2024_DEBUG_PRINTF("pace2024::branch_and_cut\tconstructing lp\n");
        construct_lp();

        // compute a first heuristic solution for an upper bound
        median_heuristic(graph, ordering).run();
        upper_bound = crossing_number_of<T, R>(graph, ordering);
    }

    // delete copy constructor, move constructor, copy assignment and move assignment
    branch_and_cut(const branch_and_cut &rhs) = delete;
    branch_and_cut(branch_and_cut &&rhs) = delete;
    branch_and_cut &operator=(const branch_and_cut &rhs) = delete;
    branch_and_cut &operator=(branch_and_cut &&rhs) = delete;

    /**
     * @brief deletes the lp solver object
     */
    ~branch_and_cut() {
        glp_delete_prob(lp);
    }

   private:
    //
    // initializer methods
    //

    /**
     * @brief constructs the lp for branch and cut:
     * - computes the crossing numbers
     * - adds the corresponding variables to the lp
     * - fixes relative positions if possible
     */
    inline void construct_lp() {
        int k = glp_add_cols(lp, static_cast<int>(n1_choose_2));
        for (T i = 0; i < n1; ++i) {
            for (T j = i + 1; j < n1; ++j) {
                add_variable_to_lp(i, j, k);
                ++k;
            }
        }
        // set constant term (shift/offset) in the objective function
        glp_set_obj_coef(lp, 0, static_cast<double>(obj_val_offset));
    }

    /**
     * @brief computes the crossing number, sets the corresponding bounds for the variable,
     * fixes it if possible, and computes lower_bound and obj_val_offsets (partially).
     *
     * @param i vertex i
     * @param j vertex j
     * @param k index k for lp
     */
    inline void add_variable_to_lp(const T i, const T j, const int k) {
        assert(static_cast<std::size_t>(k) <= n1_choose_2);
        auto [c_ij, c_ji] = crossing_numbers_of<T, R>(graph, i, j);

        // lower_bound = sum min(c_ij, c_ji)
        if (c_ij < c_ji) {
            lower_bound += c_ij;
            initial_solution[k] = 1;
        } else {
            lower_bound += c_ji;
            initial_solution[k] = 0;
        }

        // obj_val_offset = sum c_ji (since we subtract the c_ji for the coefficients)
        obj_val_offset += c_ji;

        if (c_ij == 0) {
            // fix i < j in the ordering
            glp_set_col_bnds(lp, k, GLP_FX, 1., 1.);  // ub is ignored
        } else if (c_ji == 0) {
            // fix j < i in the ordering
            glp_set_col_bnds(lp, k, GLP_FX, 0., 0.);  // ub is ignored
        } else {
            // set 0 <= x_ij <= 1
            glp_set_col_bnds(lp, k, GLP_DB, 0., 1.);
        }
        // todo: if d(i)=d(j), check if fixing is possible

        // set coefficient of added column
        double coef = static_cast<double>(c_ij) - static_cast<double>(c_ji);
        glp_set_obj_coef(lp, k, coef);
    }

    //
    // helper methods
    //

   public:
    /**
     * @brief converts an index pair i, j, i < j, to the lp column (variable)
     * index
     *
     * @param i
     * @param j
     * @return int lp column index
     */
    inline int get_variable_index(const int i, const int j) {
        assert(i < j);
        int offset = n1_choose_2 - (n1 - i) * (n1 - i - 1) / 2;
        int index = offset + j - i;
        assert(1 <= index && index <= static_cast<int>(n1_choose_2));
        return index;
    }

   private:
    /**
     * @brief given by how much the 3-cycle ieq is violated, that is,
     * # x_ij + x_jk - x_ik - 1 > 0
     * # -x_ij - x_jk + x_ik > 0
     * returns the corresponding bucket index
     *
     * @param val in (0,1]
     * @return std::size_t bucket index
     */
    inline std::size_t get_bucket(double val) {
        assert(0 < val);
        assert(val <= 1);
        const std::size_t i = static_cast<std::size_t>(val * (PACE2024_CONST_NOF_BUCKETS - 1));
        assert(i < PACE2024_CONST_NOF_BUCKETS);
        return i;
    }

    /**
     * @brief checks if a column (variable) of the lp is integral
     *
     * @param j index between 1 and n1_choose_2
     * @return true if integral
     * @return false if not
     */
    inline bool is_column_integral(int j) {
        double x = glp_get_col_prim(lp, j);
        if (x > PACE2024_CONST_EPSILON && x < 1. - PACE2024_CONST_EPSILON) {
            return false;
        }
        return true;
    }

    /**
     * @brief checks if the current solution is integral
     *
     * @return int 0, if solution is integral, index j, 1 <= j <= n1_choose_2,
     * of nonintegral column if not
     */
    int is_solution_integral() {
        for (int j = 1; j <= static_cast<int>(n1_choose_2); ++j) {
            if (!is_column_integral(j)) {
                return j;
            }
        }
        return 0;
    }

    /**
     * @brief constructs the restrictions graph into digraph and
     * then computes a topological sort of digraph for the ordering.
     * stores the ordering in ordering and the number of crossings in upper_bound.
     * requirement: lp solution integral.
     */
    void compute_ordering() {
        // resets the graph
        digraph.rollback();

        // build the constraint graph
        int k = 1;
        for (T i = 0; i < n1; ++i) {
            for (T j = i + 1; j < n1; ++j) {
                assert(is_column_integral(k));
                double x = glp_get_col_prim(lp, k);
                if (x < PACE2024_CONST_EPSILON) {
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
        double value = glp_get_obj_val(lp);
        assert(value >= 0);
        upper_bound = static_cast<R>(llround(value));
    }

    //
    // cutting plane methods
    //

    /**
     * @brief adds the forward or backward 3-cycle ieq (depending on forward)
     * for the variables with indices ij, jk and ik
     *
     * @tparam forward
     * @param ij
     * @param jk
     * @param ik
     */
    template <bool forward>
    inline void add_3cycle_ieq(const int &ij, const int &jk, const int &ik) {
        if constexpr (forward) {
            int row = glp_add_rows(lp, 1);
            glp_set_row_bnds(lp, row, GLP_UP, 0., 1.);
            const int indices[4] = {0, ij, jk, ik};
            const double coefficients[4] = {0, 1., 1., -1.};
            glp_set_mat_row(lp, row, 3, indices, coefficients);
        } else {
            int row = glp_add_rows(lp, 1);
            glp_set_row_bnds(lp, row, GLP_UP, 0., 0.);
            const int indices[4] = {0, ij, jk, ik};
            const double coefficients[4] = {0, -1., -1., 1.};
            glp_set_mat_row(lp, row, 3, indices, coefficients);
        }
    }

    /**
     * @brief given the row indices ij, jk and ik, returns
     * # x_ij + x_jk - x_ik if forward = true
     * # -x_ij - x_jk + x_ik if forward = false
     *
     * @tparam forward determines if we consider the forward or backward cycle
     * @param ij
     * @param jk
     * @param ik
     * @return double
     */
    template <bool forward>
    inline double get_3cycle_ieq_value(const int &ij, const int &jk, const int &ik) {
        double x_ij = glp_get_col_prim(lp, ij);
        double x_jk = glp_get_col_prim(lp, jk);
        double x_ik = glp_get_col_prim(lp, ik);

        if constexpr (forward) {
            return x_ij + x_jk - x_ik;
        } else {
            return -x_ij - x_jk + x_ik;
        }
    }

    /**
     * @brief given the vertices i, j, k, checks the 3-cycle constraints
     * ijk and ikj. if violated, it adds them to the lp, too.
     *
     * @param i
     * @param j
     * @param k
     * @return std::size_t number of found and added 3-cycle constraints: 0, 1 or 2
     */
    inline std::size_t check_3cycle(const int i, const int j, const int k) {
        const int ij = get_variable_index(i, j),
                  jk = get_variable_index(j, k),
                  ik = get_variable_index(i, k);

        std::size_t nof_new_cycle_constraints = 0;

        // cycle ijk
        if (get_3cycle_ieq_value<true>(ij, jk, ik) > 1) {
            ++nof_new_cycle_constraints;
            add_3cycle_ieq<true>(ij, jk, ik);
        }

        // cycle ikj
        if (get_3cycle_ieq_value<false>(ij, jk, ik) > 0) {
            ++nof_new_cycle_constraints;
            add_3cycle_ieq<false>(ij, jk, ik);
        }

        return nof_new_cycle_constraints;
    }

    /**
     * @brief tries to find violated 3-cycle constraints
     * and adds the first <= PACE2024_CONST_NOF_CYCLE_CONSTRAINTS to the lp
     *
     * @return true 3-cycle found
     * @return false 3-cycle not found
     */
    bool check_3cycles() {
        std::size_t nof_cycle_constraints = 0;

        for (T i = 0; i < n1; ++i) {
            for (T j = i + 1; j < n1; ++j) {
                for (T k = j + 1; k < n1; ++k) {
                    nof_cycle_constraints += check_3cycle(i, j, k);

                    if (nof_cycle_constraints >
                        PACE2024_CONST_NOF_CYCLE_CONSTRAINTS) {
                        return true;
                    }
                }
            }
        }

        return nof_cycle_constraints > 0;
    }

    bool check_3cycle_and_bucket_it(const int i, const int j, const int k) {
        const int ij = get_variable_index(i, j),
                  jk = get_variable_index(j, k),
                  ik = get_variable_index(i, k);
        return false;
    }

    bool check_3cycles_exhaustively() {
        return false;
    }

    /**
     * @brief search for violated constraints
     * and add these to the lp.
     *
     * @return true some where found
     * @return false otherwise
     */
    bool try_to_generate_cutting_planes() {
        PACE2024_DEBUG_PRINTF("pace2024::branch_and_cut\t\tstart checking cycle constraints\n");
        bool success = check_3cycles();
        PACE2024_DEBUG_PRINTF("pace2024::branch_and_cut\t\tend checking cycle constraints (%s)\n", success);
        return success;
    }

    //
    // variable fixing methods
    //

    /**
     * @brief fixes variables according to the reduced cost condition
     */
    void perform_permanent_fixing() {
        for (int j = 1; j <= static_cast<int>(n1_choose_2); ++j) {
            // int col_stat = glp_get_col_stat(lp, j);
            // PACE2024_DEBUG_PRINTF("col_stat=%s\n", col_stat);
            // if (col_stat == GLP_BS) continue;

            double reduced_cost = glp_get_col_dual(lp, j);

            if (initial_solution[j] == 0 &&
                static_cast<double>(lower_bound) + reduced_cost >= static_cast<double>(upper_bound)) {
                glp_set_col_bnds(lp, j, GLP_FX, 0., 0.);
                PACE2024_DEBUG_PRINTF("pace2024::branch_and_cut\t\tpermanent fixing of %d\n", j);
            }

            if (initial_solution[j] == 1 &&
                static_cast<double>(lower_bound) - reduced_cost >= static_cast<double>(upper_bound)) {
                glp_set_col_bnds(lp, j, GLP_FX, 1., 1.);
                PACE2024_DEBUG_PRINTF("pace2024::branch_and_cut\t\tpermanent fixing of %d\n", j);
            }
        }
    }

    //
    // branch and bound
    //

    /**
     * @brief if there's a variable on the stack,
     * the function pops it, and we fix the opposite value.
     *
     * @return true if stack is empty
     * @return false otherwise
     */
    bool backtrack_fix_column() {
        if (stack.empty()) {
            return true;
        }

        // fix the opposite value
        int j = stack.top();
        stack.pop();
        double x = glp_get_col_prim(lp, j);
        glp_set_col_bnds(lp, j, GLP_FX, 1. - x, 0.);  // ub is ignored

        return false;
    }

    /**
     * @brief temporarily fixes column j
     *
     * @param j
     */
    void fix_column(int j) {
        // todo: more sophisticated fixing
        // todo: fix implications, too
        glp_set_col_bnds(lp, j, GLP_FX, 0., 0.);
        stack.emplace(j);
    }

    /**
     * @brief given the optimal value of the current lp,
     * decides if we generate cutting planes,
     * or if we branch,
     * or if we backtrack.
     *
     * @param value optimal value of current lp
     * @return true optimal solution found
     * @return false otherwise
     */
    inline bool branch_n_cut(double value) {
        if (static_cast<double>(upper_bound) <= value) {
            // try to backtrack
            // if this is not possible (meaning that the stack is empty),
            // then we found an optimal solution
            return backtrack_fix_column();
        } else {
            // try to generate cutting planes
            bool successul = try_to_generate_cutting_planes();

            if (!successul) {
                int j = is_solution_integral();
                if (j == 0) {
                    // the solution is integral; we found a better solution
                    compute_ordering();
                    perform_permanent_fixing();
                    return backtrack_fix_column();
                } else {
                    // todo: better selection of 3-cycle ieqs to add
                    // todo: improve solution with heuristic
                    // todo: permanent fixing
                    // todo: transitivity
                    // todo: delete "loose" inequalities
                    // todo: better selection
                    // branch by fixing a column
                    fix_column(j);
                }
            }
        }
        return false;
    }

   public:
    /**
     * @brief solves the given instance exactly with a
     * branch and cut approach.
     *
     * @param do_print
     */
    void solve(bool do_print = true) {
        // solve the lp
        glp_simplex(lp, &params);
        int status = glp_get_status(lp);
        (void)status;
        assert(status == GLP_OPT);

        // perform initial permanent fixing since we already have a heuristic solution
        perform_permanent_fixing();

        // get value of current optimal solution and call branch_n_cut with it
        double value = glp_get_obj_val(lp);
        assert(static_cast<R>(llround(value)) == lower_bound);
        bool optimum = branch_n_cut(value);

        if (!optimum) {
            for (std::size_t iteration = 0;; ++iteration) {
                // solve the lp
                PACE2024_DEBUG_PRINTF("pace2024::branch_and_cut\tstart simplex\n");
                glp_simplex(lp, &params);
                status = glp_get_status(lp);
                assert(status == GLP_OPT);
                PACE2024_DEBUG_PRINTF("pace2024::branch_and_cut\tend simplex\n");

                // get value of current optimal solution and call branch_n_cut with it
                value = glp_get_obj_val(lp);
                optimum = branch_n_cut(value);
                if (optimum) break;
            }
        }

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