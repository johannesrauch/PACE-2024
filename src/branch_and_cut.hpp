#ifndef PACE2024_BRANCH_AND_CUT_HPP
#define PACE2024_BRANCH_AND_CUT_HPP

#ifndef PACE2024_CONST_NOF_CYCLE_CONSTRAINTS
#define PACE2024_CONST_NOF_CYCLE_CONSTRAINTS 1024
#endif

#ifndef PACE2024_CONST_NOF_BUCKETS
#define PACE2024_CONST_NOF_BUCKETS 10
#endif

#include <glpk.h>
#include <math.h>

#include <algorithm>
#include <cfloat>
#include <functional>
#include <stack>
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

    struct bucket_entry {
        const int ij, jk, ik;
        const bool ub;

        explicit bucket_entry(const int ij, const int jk, const int ik, const bool ub)
            : ij(ij), jk(jk), ik(ik), ub(ub) {}
    };

    /**
     * @brief buckets for bucket sorting violated 3-cycle ieqs
     */
    std::vector<std::vector<bucket_entry>> buckets{PACE2024_CONST_NOF_BUCKETS};

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
        PACE2024_DEBUG_PRINTF("start constructing lp\n");
        construct_lp();
        PACE2024_DEBUG_PRINTF("end constructing lp\n");

        // compute a first heuristic solution for an upper bound
        PACE2024_DEBUG_PRINTF("start heuristic\n");
        probabilistic_median_heuristic<T, R>(graph, ordering).run();
        PACE2024_DEBUG_PRINTF("end heuristic\n");
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
            glp_set_col_bnds(lp, k, GLP_FX, 1., 0.);  // ub is ignored
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
    // solution methods
    //

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
        double value = glp_get_obj_val(lp);
        assert(value >= 0);
        upper_bound = static_cast<R>(llround(value));
    }

    //
    // helper methods
    //

   private:
    /**
     * @brief erases contents of each bucket in buckets
     */
    inline void clear_buckets() {
        for (auto &bucket : buckets) {
            bucket.clear();
        }
    }

   public:
    /**
     * @brief converts an index pair i, j, i < j, to the lp column (variable)
     * index
     *
     * @param i
     * @param j
     * @return int lp column index
     */
    inline int get_variable_index(const int &i, const int &j) {
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
    inline std::size_t get_bucket(const double &val) {
        assert(0 < val);
        assert(val <= 1 + 3 * params.tol_bnd);
        const std::size_t i = static_cast<std::size_t>(val * (PACE2024_CONST_NOF_BUCKETS - 1));
        assert(i < PACE2024_CONST_NOF_BUCKETS);
        return i;
    }

    /**
     * @brief returns x_ij + x_jk - x_ik of the lp
     */
    inline double get_3cycle_ieq_value(const int &ij, const int &jk, const int &ik) {
        const double x_ij = glp_get_col_prim(lp, ij);
        const double x_jk = glp_get_col_prim(lp, jk);
        const double x_ik = glp_get_col_prim(lp, ik);
        return x_ij + x_jk - x_ik;
    }

    inline bool has_row_lb_slack(const int &i) {
        const double x = glp_get_row_prim(lp, i);
        return has_row_lb(i) && x > params.tol_bnd;
    }

    inline bool has_row_ub_slack(const int &i) {
        const double x = glp_get_row_prim(lp, i);
        return has_row_ub(i) && x < 1. - params.tol_bnd;
    }

    /**
     * @brief returns lb != -DBL_MAX
     */
    inline bool has_row_lb(const int &i) {
        const double lb = glp_get_row_lb(lp, i);
        return lb != -DBL_MAX;
    }

    /**
     * @brief returns ub != DBL_MAX
     */
    inline bool has_row_ub(const int &i) {
        const double ub = glp_get_row_ub(lp, i);
        return ub != DBL_MAX;
    }

    /**
     * @brief returns x < -3 * 1e-7
     * (the factor 3 comes from the sum x_ij + x_jk - x_ik, where each factor has a tolerance of 1e-7)
     */
    inline bool is_3cycle_lb_violated(const double &x) {
        return x < -3 * params.tol_bnd;
    }

    /**
     * @brief returns 1 + 3 * 1e-7
     * (the factor 3 comes from the sum x_ij + x_jk - x_ik, where each factor has a tolerance of 1e-7)
     */
    inline bool is_3cycle_ub_violated(const double &x) {
        return x > 1. + 3 * params.tol_bnd;
    }

    /**
     * @brief returns if the last bucket is full
     *
     * @return true if last bucket has >= PACE2024_CONST_NOF_CYCLE_CONSTRAINTS elements
     * @return false otherwise
     */
    inline bool is_last_bucket_full() {
        return (*buckets.rbegin()).size() >= PACE2024_CONST_NOF_CYCLE_CONSTRAINTS;
    }

    /**
     * @brief checks if a column (variable) of the lp is integral
     *
     * @param j index between 1 and n1_choose_2
     * @return true if integral
     * @return false if not
     */
    inline bool is_column_integral(const int &j) {
        double x = glp_get_col_prim(lp, j);
        if (x > params.tol_bnd && x < 1. - params.tol_bnd) {
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

    //
    // cutting plane methods
    //

    inline int add_3cycle_row(const int &ij, const int &jk, const int &ik) {
        const int row = glp_add_rows(lp, 1);
        const int indices[4] = {0, ij, jk, ik};
        const double coefficients[4] = {0, 1., 1., -1.};
        glp_set_mat_row(lp, row, 3, indices, coefficients);
        return row;
    }

    inline void add_3cycle_ieq_lb(const int &ij, const int &jk, const int &ik) {
        const int row = add_3cycle_row(ij, jk, ik);
        glp_set_row_bnds(lp, row, GLP_LO, 0., 1.);
    }

    inline void add_3cycle_ieq_ub(const int &ij, const int &jk, const int &ik) {
        const int row = add_3cycle_row(ij, jk, ik);
        glp_set_row_bnds(lp, row, GLP_UP, 0., 1.);
    }

    /**
     * @brief expects the violated 3-cycle ieqs in buckets.
     * adds the most violated <= PACE2024_CONST_NOF_CYCLE_CONSTRAINTS to the lp.
     *
     * @return std::size_t number of new rows
     */
    inline std::size_t add_3cycle_ieq_bnds() {
        std::size_t nof_new_rows{0};

        for (auto r_it = buckets.rbegin(); r_it != buckets.rend(); ++r_it) {
            for (const auto &[ij, jk, ik, ub] : *r_it) {
                ++nof_new_rows;

                if (ub) {
                    add_3cycle_ieq_ub(ij, jk, ik);
                } else {
                    add_3cycle_ieq_lb(ij, jk, ik);
                }

                if (nof_new_rows >= PACE2024_CONST_NOF_CYCLE_CONSTRAINTS) {
                    return nof_new_rows;
                }
            }
        }

        return nof_new_rows;
    }

    inline void check_3cycle_ieq_lb(const int &ij, const int &jk, const int &ik, const double &x) {
        if (is_3cycle_lb_violated(x)) {
            buckets[get_bucket(-x)].emplace_back(ij, jk, ik, false);
        }
    }

    inline void check_3cycle_ieq_ub(const int &ij, const int &jk, const int &ik, const double &x) {
        if (is_3cycle_ub_violated(x)) {
            buckets[get_bucket(x - 1)].emplace_back(ij, jk, ik, true);
        }
    }

    inline void check_3cycle(const T &i, const T &j, const T &k) {
        const int ij = get_variable_index(i, j);
        const int jk = get_variable_index(j, k);
        const int ik = get_variable_index(i, k);
        const double x = get_3cycle_ieq_value(ij, jk, ik);
        check_3cycle_ieq_lb(ij, jk, ik, x);
        check_3cycle_ieq_ub(ij, jk, ik, x);
    }

    bool check_3cycles() {
        PACE2024_DEBUG_PRINTF("\tstart check_3cycles\n");
        
        clear_buckets();
        for (T i = 0; i < n1; ++i) {
            for (T j = i + 1; j < n1; ++j) {
                for (T k = j + 1; k < n1; ++k) {
                    check_3cycle(i, j, k);
                }
            }
        }

        const std::size_t nof_new_rows = add_3cycle_ieq_bnds();
        PACE2024_DEBUG_PRINTF("\tend check_3cycles, number of new rows=%lld\n", nof_new_rows);
        return nof_new_rows > 0;
    }

    /**
     * @brief search for violated constraints
     * and add these to the lp.
     *
     * @return true some where found
     * @return false otherwise
     */
    bool cut() {
        bool success = check_3cycles();
        return success;
    }

    //
    // ieq deletion methods
    //

    /**
     * @brief delete ieqs with postive slack from the lp (if stack.size() > 0)
     */
    void remove_positive_slack_ieqs() {
        if (stack.size() > 0) return;
        PACE2024_DEBUG_PRINTF("\tstart remove_positive_slack_ieqs\n");

        std::vector<int> rows_to_remove;
        rows_to_remove.emplace_back(0);

        // gather rows to delete in respective vector or just remove lb/ub's from row
        const int nof_rows = glp_get_num_rows(lp);
        for (int i = 1; i <= nof_rows; ++i) {
            if (has_row_lb(i)) {
                if (has_row_lb_slack(i)) {
                    rows_to_remove.emplace_back(i);
                }
            } else {
                assert(has_row_ub(i));
                if (has_row_ub_slack(i)) {
                    rows_to_remove.emplace_back(i);
                }
            }
        }

        const int nof_removed_rows = static_cast<int>(rows_to_remove.size()) - 1;
        if (nof_removed_rows > 0) {
            PACE2024_DEBUG_PRINTF("\tend remove_positive_slack_ieqs, number of removed rows=%lld\n", nof_removed_rows);
            glp_del_rows(lp, nof_removed_rows, &rows_to_remove[0]);
        }
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
                PACE2024_DEBUG_PRINTF("\tpermanent fixing of %d\n", j);
            }

            if (initial_solution[j] == 1 &&
                static_cast<double>(lower_bound) - reduced_cost >= static_cast<double>(upper_bound)) {
                glp_set_col_bnds(lp, j, GLP_FX, 1., 0.);
                PACE2024_DEBUG_PRINTF("\tpermanent fixing of %d\n", j);
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
    bool backtrack() {
        if (stack.empty()) {
            return true;
        }

        // fix the opposite value
        const int j = stack.top();
        stack.pop();
        const double fix = glp_get_col_lb(lp, j);
        if (fix > 0.5) {
            PACE2024_DEBUG_PRINTF("\tfix_column %d to 0\n", j);
            glp_set_col_bnds(lp, j, GLP_FX, 0., 0.);  // ub is ignored
        } else {
            PACE2024_DEBUG_PRINTF("\tfix_column %d to 1\n", j);
            glp_set_col_bnds(lp, j, GLP_FX, 1., 0.);  // ub is ignored
        }

        return false;
    }

    /**
     * @brief temporarily fixes column j
     *
     * @param j
     */
    inline void branch(const int &j) {
        // todo: more sophisticated fixing
        // todo: fix implications, too
        PACE2024_DEBUG_PRINTF("\tfix_column %d to 0\n", j);
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
    inline bool branch_n_cut() {
        double value = glp_get_obj_val(lp);

        if (static_cast<double>(upper_bound) <= value) {
            // try to backtrack
            // if this is not possible (meaning that the stack is empty),
            // then we found an optimal solution
            return backtrack();
        } else {
            // try to generate cutting planes
            bool successul = cut();

            if (!successul) {
                int j = is_solution_integral();
                if (j == 0) {
                    // the solution is integral; we found a better solution
                    compute_ordering();
                    perform_permanent_fixing();
                    return backtrack();
                } else {
                    // todo: improve solution with heuristic
                    // todo: transitivity
                    // todo: better selection
                    // branch by fixing a column
                    branch(j);
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
        bool is_optimum = branch_n_cut();
        assert(get_rounded_obj_value() == lower_bound);

        // driver loop
        for (std::size_t iteration = 0; !is_optimum; ++iteration) {
            PACE2024_DEBUG_PRINTF("depth=%d, upper_bound=%d, nof_rows=%d\n",
                                  stack.size(),
                                  upper_bound,
                                  glp_get_num_rows(lp));

            // solve the lp
            PACE2024_DEBUG_PRINTF("start glp_simplex\n", stack.size());
            glp_simplex(lp, &params);
            status = glp_get_status(lp);
            (void) status;
            // assert(status == GLP_OPT);
            PACE2024_DEBUG_PRINTF("end glp_simplex, status=%d, objective value=%f\n", status, glp_get_obj_val(lp));

            // delete untight ieqs
            remove_positive_slack_ieqs();

            // branch or cut
            is_optimum = branch_n_cut();
        }

        // output
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

   private:
    /**
     * @brief retrieve objective value from lp and round it
     *
     * @return R rounded objective value of lp
     */
    inline R get_rounded_obj_value() {
        return static_cast<R>(llround(glp_get_obj_val(lp)));
    }
};

};  // namespace pace2024

#endif