#ifndef PACE2024_BRANCH_AND_CUT_HPP
#define PACE2024_BRANCH_AND_CUT_HPP

#ifndef PACE2024_CONST_NOF_CYCLE_CONSTRAINTS
#define PACE2024_CONST_NOF_CYCLE_CONSTRAINTS 8
#endif

// 1e-7 is the default tolerance in glpk for feasibility
#ifndef PACE2024_CONST_EPSILON
#define PACE2024_CONST_EPSILON 1e-7
#endif

#include <glpk.h>

#include <stack>
#include <math.h>

#include "bipartite_graph.hpp"
#include "crossings.hpp"
#include "instance.hpp"
#include "matrix.hpp"
#include "median_heuristic.hpp"
#include "topological_sort.hpp"

namespace pace2024 {

template <typename T>
class branch_and_cut {
   private:
    // instance
    const general_bipartite_graph<T> graph;
    const int n1, n1_choose_2;
    const folded_square_matrix<T> cr_matrix;

    // lp
    glp_prob* const lp;
    glp_smcp params;

    // solution
    uint64_t lower_bound;
    uint64_t upper_bound;  // current best solution
    std::vector<T> ordering;

    // internal variables
    std::stack<int> stack;

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
    branch_and_cut(const general_instance<T>& instance)
        : graph(instance),
          n1(graph.get_n1()),
          n1_choose_2(n1 * (n1 - 1) / 2),
          cr_matrix(graph),
          lp(glp_create_prob()),
          lower_bound(0),
          upper_bound(0),
          ordering(n1) {
        glp_init_smcp(&params);
        glp_set_obj_dir(lp, GLP_MIN);

        // add n1(n1-1)/2 variables
        glp_add_cols(lp, n1_choose_2);
        int k = 1;
        for (int i = 0; i < n1; ++i) {
            for (int j = i + 1; j < n1; ++j) {
                // just a test
                assert(k == get_variable_index(i, j));

                T c_ij = cr_matrix(i, j), c_ji = cr_matrix(j, i),
                  c = c_ij - c_ji;
                lower_bound += c_ij < c_ji ? c_ij : c_ji;  // min(c_ij, c_ji)

                if (c_ij == 0) {
                    glp_set_col_bnds(lp, k, GLP_FX, 1, 1);  // fix i < j
                } else if (c_ji == 0) {
                    glp_set_col_bnds(lp, k, GLP_FX, 0, 0);  // fix j < i
                } else {
                    // set 0 <= x_ij <= 1
                    glp_set_col_bnds(lp, k, GLP_DB, 0, 1);
                }

                glp_set_obj_coef(lp, k, (double)c);
                ++k;
            }
        }

        median_heuristic(graph, ordering);
        upper_bound = compute_crossings(cr_matrix, ordering);
    }

    ~branch_and_cut() { glp_delete_prob(lp); }

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
        assert(index <= n1_choose_2);
        return index;
    }

   private:
    /**
     * @brief given indices i, j, k, checks the 3-cycle constraints ijk and ikj
     * if violated, it adds them to the lp, too
     *
     * @param i
     * @param j
     * @param k
     * @return T number of found and added 3-cycle constraints: 0, 1
     * or 2
     */
    inline T check_cycle_constraint(const int i, const int j, const int k) {
        const int index_ij = get_variable_index(i, j),
                  index_jk = get_variable_index(j, k),
                  index_ik = get_variable_index(i, k);

        double x_ij = glp_get_col_prim(lp, index_ij);
        double x_jk = glp_get_col_prim(lp, index_jk);
        double x_ik = glp_get_col_prim(lp, index_ik);

        std::size_t nof_new_cycle_constraints = 0;
        // cycle ijk
        if (x_ij + x_jk - x_ik > 1) {
            ++nof_new_cycle_constraints;
            int row = glp_add_rows(lp, 1);
            glp_set_row_bnds(lp, row, GLP_UP, 0., 1.);
            const int indices[4] = {0, index_ij, index_jk, index_ik};
            const double coefficients[4] = {0, 1., 1., -1.};
            glp_set_mat_row(lp, row, 3, indices, coefficients);
        }

        // cycle ikj
        if (x_ik - x_ij - x_jk > 0) {
            ++nof_new_cycle_constraints;
            int row = glp_add_rows(lp, 1);
            glp_set_row_bnds(lp, row, GLP_UP, 0., 0.);
            const int indices[4] = {0, index_ij, index_jk, index_ik};
            const double coefficients[4] = {0, -1., -1., 1.};
            glp_set_mat_row(lp, row, 3, indices, coefficients);
        }

        return nof_new_cycle_constraints;
    }

    /**
     * @brief tries to find violated 3-cycle constraints
     * and adds them to the lp
     *
     * @return true 3-cycle found
     * @return false 3-cycle not found
     */
    bool check_cycle_constraints() {
        T nof_cycle_constraints = 0;

        for (int i = 0; i < n1; ++i) {
            for (int j = 0; j < n1; ++j) {
                for (int k = 0; k < n1; ++k) {
                    nof_cycle_constraints += check_cycle_constraint(i, j, k);

                    if (nof_cycle_constraints >
                        PACE2024_CONST_NOF_CYCLE_CONSTRAINTS) {
                        return true;
                    }
                }
            }
        }

        return nof_cycle_constraints > 0;
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
        for (int j = 1; j <= n1_choose_2; ++j) {
            if (!is_column_integral(j)) {
                return j;
            }
        }
        return 0;
    }

    bool try_to_generate_cutting_planes() { return check_cycle_constraints(); }

    void compute_ordering() {
        // build the constraint graph
        general_graph<T> graph(n1_choose_2);
        int k = 1;
        for (T i = 0; i < n1; ++i) {
            for (T j = 0; j < n1; ++j) {
                assert(is_column_integral(k));
                double x = glp_get_col_prim(lp, k);
                if (x < PACE2024_CONST_EPSILON) {
                    graph.add_edge(j, i);
                } else {
                    graph.add_edge(i, j);
                }
                ++k;
            }
        }

        // the ordering computed by the lp is the topological sort
        assert(topological_sort(graph, ordering));
        upper_bound = llround(glp_get_obj_val(lp));
    }

    /**
     * @brief if there's a variable on the stack,
     * the function pops it, and we fix the opposite value
     *
     * @return true if stack is empty
     * @return false if not
     */
    bool backtrack() {
        if (stack.empty()) {
            return true;
        }

        // fix the opposite value
        int j = stack.top();
        stack.pop();
        double x = glp_get_col_prim(lp, j);
        glp_set_col_bnds(lp, j, GLP_FX, 1. - x, 1. - x);

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

   public:
    void solve() {
        for (;;) {
            glp_simplex(lp, &params);
            int status = glp_get_status(lp);
            assert(status == GLP_OPT);

            double value = glp_get_obj_val(lp);
            if ((double)upper_bound <= value) {
                if (backtrack()) return;
            } else {
                if (!try_to_generate_cutting_planes()) {
                    int j = is_solution_integral();
                    if (j == 0) {
                        // we found a better solution
                        compute_ordering();
                        // todo: permanent fixing
                    } else {
                        // todo: improve solution with heuristic
                        // todo: permanent fixing
                        // todo: better selection
                        fix_column(j);
                    }
                }
            }
        }
    }
};

};  // namespace pace2024

#endif