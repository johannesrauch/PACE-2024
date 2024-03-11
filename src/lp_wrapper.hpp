#ifndef PACE2024_LP_WRAPPER_HPP
#define PACE2024_LP_WRAPPER_HPP

#include <glpk.h>
#include <highs/Highs.h>

#include <cassert>
#include <coin/ClpModel.hpp>
#include <coin/ClpSimplex.hpp>
#include <coin/CoinBuild.hpp>

#include "bipartite_graph.hpp"
#include "crossing_number.hpp"
#include "debug_printf.hpp"

#ifndef PACE2024_CONST_NOF_CYCLE_CONSTRAINTS
#define PACE2024_CONST_NOF_CYCLE_CONSTRAINTS 512
#endif

#ifndef PACE2024_CONST_NOF_BUCKETS
#define PACE2024_CONST_NOF_BUCKETS 10
#endif

#ifndef PACE2024_CONST_FEASIBILITY_TOLERANCE
#define PACE2024_CONST_FEASIBILITY_TOLERANCE 1e-7
#endif

#ifndef PACE2024_CONST_INTEGER_TOLERANCE
#define PACE2024_CONST_INTEGER_TOLERANCE 1e-6
#endif

namespace pace2024 {

/**
 * @brief generic wrapper class for an lp solver used by branch_and_cut
 */
class lp_wrapper {
   protected:
    //
    // bucket attributes and methods
    //

    /// @brief struct used for storing violated 3-cycle inequalities
    struct bucket_entry {
        /// @brief column indices
        const int ij, jk, ik;
        /// @brief true is 3-cycle ijk, false is 3-cycle ikj
        const bool ub;

        explicit bucket_entry(const int ij, const int jk, const int ik, const bool ub)
            : ij(ij), jk(jk), ik(ik), ub(ub) {}
    };

    /// @brief buckets for bucket sorting violated 3-cycle inequalities
    std::vector<std::vector<bucket_entry>> buckets{PACE2024_CONST_NOF_BUCKETS};

    /// @brief erases contents of each bucket in buckets
    inline void clear_buckets() {
        for (auto &bucket : buckets) {
            bucket.clear();
        }
    }

    /**
     * @brief given by how much the 3-cycle inequality is violated,
     * returns the corresponding bucket
     *
     * @param val in (0,1]
     * @return std::vector<bucket_entry>& the corresponding bucket
     */
    inline std::vector<bucket_entry> &get_bucket(const double &val) {
        assert(0 < val);
        assert(val <= 1);
        const int i = static_cast<int>(val * PACE2024_CONST_NOF_BUCKETS - 1);
        assert(0 <= i);
        assert(i < PACE2024_CONST_NOF_BUCKETS);
        return buckets[i];
    }

    /// @brief returns sum of the number of elements in each bucket of `buckets`
    inline int get_nof_bucket_entries() {
        int n{0};
        for (const auto &bucket : buckets) {
            n += bucket.size();
        }
        return n;
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

    //
    // lp related attributes
    //

    /// @brief best lower bound
    int lower_bound{0};

    /// @brief best upper bound
    int upper_bound{INT_MAX};

    //
    // instance related attributes
    //

    /// @brief number of vertices in free layer
    const int n1;

    /// @brief `n1_choose_2` = `n1` * (`n1` - 1) / 2
    const int n1_choose_2;

   public:
    lp_wrapper(const int n1) : n1(n1), n1_choose_2(n1 * (n1 - 1) / 2) {}
    virtual ~lp_wrapper() {}
    virtual bool cut() = 0;
    virtual void fix_column(const int j, const double fix_to) = 0;
    virtual void fix_columns(const int new_upper_bound) = 0;
    virtual double get_column_value(const int j) = 0;
    virtual int get_nof_cols() = 0;
    virtual int get_nof_rows() = 0;
    virtual double get_objective_value() = 0;
    virtual int64_t get_rounded_objective_value() = 0;
    virtual int is_integral() = 0;
    virtual bool is_optimal() = 0;
    virtual void solve(bool delete_rows_after) = 0;
    virtual void unfix_column(const int j) = 0;
};

//
// wrapper class for glpk solver
//

class glpk_wrapper : public lp_wrapper {
    /// @brief const pointer to (mutable) lp solver
    glp_prob *const lp;
    /// @brief parameters for glp_simplex
    glp_smcp params;

   public:
    template <typename T>
    glpk_wrapper(const general_bipartite_graph<T> &graph,
                 const int msg_level = GLP_MSG_OFF)
        : lp_wrapper(static_cast<int>(graph.get_n1())),
          lp(glp_create_prob()) {
        initialize_parameters(msg_level);
        add_columns(graph);
        glp_set_obj_dir(lp, GLP_MIN);
    }

    glpk_wrapper(const glpk_wrapper &rhs) = delete;
    glpk_wrapper(glpk_wrapper &&rhs) = delete;
    glpk_wrapper &operator=(const glpk_wrapper &rhs) = delete;
    glpk_wrapper &operator=(glpk_wrapper &&rhs) = delete;

    ~glpk_wrapper() {
        glp_delete_prob(lp);
    }

    /**
     * @brief search for violated constraints and add these to the lp
     *
     * @return true if successful
     * @return false otherwise
     */
    virtual bool cut() {
        bool success = check_3cycles();
        // todo: k-fence
        return success;
    }

    /**
     * @brief fixes column j to fix_to
     *
     * @param j column index
     * @param fix_to
     */
    virtual void fix_column(const int j, const double fix_to) {
        assert(0 < j);
        assert(j <= glp_get_num_cols(lp));
        PACE2024_DEBUG_PRINTF("fixed variable %5d to %1.0f\n", j, fix_to);
        glp_set_col_bnds(lp, j, GLP_FX, fix_to, 0.);  // ub is ignored
    }

    /**
     * @brief fixes columns based on the following:
     * if lower_bound = sum min(c_ij, c_ji), where c_ij and c_ji are crossing numbers,
     * and coeff = c_ij - c_ji is the objective coefficient of variable x_ij,
     * and abs(c_ij - c_ji) >= upper_bound - lower_bound,
     * then we are able to fix x_ij to 0 if coeff > 0 and to 1 otherwise
     *
     * @param new_upper_bound new upper bound on the optimal value
     */
    virtual void fix_columns(const int new_upper_bound) {
        assert(new_upper_bound < upper_bound);
        upper_bound = new_upper_bound;
        for (int j = 1; j <= n1_choose_2; ++j) {
            const int diff = upper_bound - lower_bound;
            const double coeff = glp_get_obj_coef(lp, j);

            if (std::abs(coeff) >= static_cast<double>(diff) && coeff != 0.) {
                fix_column(j, coeff > 0 ? 0. : 1.);
                PACE2024_DEBUG_PRINTF("(optimality condition)\n");
            }
        }
    }

    /// @brief returns value of column j
    virtual double get_column_value(const int j) {
        return glp_get_col_prim(lp, j);
    }

    /// @brief returns the number of cols
    virtual int get_nof_cols() {
        return glp_get_num_cols(lp);
    }

    /// @brief returns the number of rows
    virtual int get_nof_rows() {
        return glp_get_num_rows(lp);
    }

    /// @brief returns the objective value of the lp
    virtual double get_objective_value() {
        return glp_get_obj_val(lp);
    }

    /// @brief returns the objective value of the lp rounded to the next integer
    virtual int64_t get_rounded_objective_value() {
        const double value = get_objective_value();
        assert(value >= 0);
        return llround(value);
    }

    /**
     * @brief returns a value based on the integrality of the current solution
     *
     * @return int 0, if solution is integral
     * @return int j, 1 <= j <= n1_choose_2, of nonintegral column otherwise
     */
    int is_integral() {
        for (int j = 1; j <= n1_choose_2; ++j) {
            if (!is_column_integral(j)) {
                return j;
            }
        }
        return 0;
    }

    /// @brief returns if an optimal feasible solution has been found
    virtual bool is_optimal() {
        return glp_get_status(lp) == GLP_OPT;
    }

    virtual void solve(bool delete_rows_after) {
        PACE2024_DEBUG_PRINTF("start glp_simplex\n");
        glp_simplex(lp, &params);
        // and since we had problems with integrality tolerances we do an exact solve, too
        glp_exact(lp, &params);
        PACE2024_DEBUG_PRINTF("end   glp_simplex, status=%d, objective value=%f\n",
                              glp_get_status(lp),
                              glp_get_obj_val(lp));

        if (delete_rows_after) {
            delete_positive_slack_rows();
        }
    }

    /// @brief resets bounds of column j to 0 <= . <= 1
    virtual void unfix_column(const int j) {
        glp_set_col_bnds(lp, j, GLP_DB, 0., 1.);
    }

   private:
    //
    // initialization methods
    //

    /**
     * @brief adds all columns of the ilp formulation of one-sided crossing minimization
     * and computes their objective coefficients (that is, the crossing numbers)
     *
     * @tparam T vertex param
     * @param graph the input graph
     */
    template <typename T>
    inline void add_columns(const general_bipartite_graph<T> &graph) {
        int obj_val_offset{0};
        int k = glp_add_cols(lp, n1_choose_2);

        for (T i = 0; i < static_cast<T>(n1); ++i) {
            for (T j = i + 1; j < static_cast<T>(n1); ++j) {
                obj_val_offset += add_column(graph, i, j, k);
                ++k;
            }
        }
        // set constant term (shift/offset) in the objective function
        glp_set_obj_coef(lp, 0, static_cast<double>(obj_val_offset));

        assert(n1_choose_2 == glp_get_num_cols(lp));
    }

    /**
     * @brief computes the crossing numbers c_ij and c_ji,
     * sets the corresponding bounds for the column (perhaps even fixes them),
     *
     * @tparam T vertex type
     * @param graph the input graph
     * @param i vertex of free layer
     * @param j vertex of free layer
     * @param k column index
     * @return crossing number c_ji (for the objective offset)
     */
    template <typename T>
    inline int add_column(const general_bipartite_graph<T> &graph,
                          const T &i,
                          const T &j,
                          const int &k) {
        assert(k <= n1_choose_2);
        const auto [c_ij, c_ji] = crossing_numbers_of<T, int>(graph, i, j);

        // lower_bound = sum min(c_ij, c_ji)
        if (c_ij < c_ji) {
            lower_bound += c_ij;
        } else {
            lower_bound += c_ji;
        }

        if (c_ij == 0 && c_ji != 0) {
            // fix i < j in the ordering
            fix_column(k, 1.);
        } else if (c_ji == 0 && c_ij != 0) {
            // fix j < i in the ordering
            fix_column(k, 0.);
        } else {
            // set 0 <= x_ij <= 1
            glp_set_col_bnds(lp, k, GLP_DB, 0., 1.);
        }

        // todo: if d(i)=d(j), check if fixing is possible

        // set coefficient of added column
        if (c_ij != c_ji) {
            const double coeff = static_cast<double>(c_ij) - static_cast<double>(c_ji);
            glp_set_obj_coef(lp, k, coeff);
        }

        return c_ji;
    }

    inline void initialize_parameters(const int &msg_level) {
        glp_init_smcp(&params);
        params.msg_lev = msg_level;
        params.meth = GLP_DUALP;
    }

    //
    // row addition methods
    //

    inline int add_3cycle_row(const int &ij, const int &jk, const int &ik) {
        const int row = glp_add_rows(lp, 1);
        const int indices[4] = {0, ij, ik, jk};
        const double coefficients[4] = {0, 1., -1., 1.};
        glp_set_mat_row(lp, row, 3, indices, coefficients);
        return row;
    }

    inline void add_3cycle_row_with_lb(const int &ij, const int &jk, const int &ik) {
        const int row = add_3cycle_row(ij, jk, ik);
        glp_set_row_bnds(lp, row, GLP_LO, 0., 1.);
    }

    inline void add_3cycle_row_with_ub(const int &ij, const int &jk, const int &ik) {
        const int row = add_3cycle_row(ij, jk, ik);
        glp_set_row_bnds(lp, row, GLP_UP, 0., 1.);
    }

    /**
     * @brief expects the violated 3-cycle ieqs in buckets.
     * adds the most violated <= PACE2024_CONST_NOF_CYCLE_CONSTRAINTS to the lp.
     *
     * @return int number of new rows
     */
    inline int add_3cycle_rows() {
        int nof_new_rows{0};

        for (auto r_it = buckets.rbegin(); r_it != buckets.rend(); ++r_it) {
            for (const auto &[ij, jk, ik, ub] : *r_it) {
                ++nof_new_rows;

                if (ub) {
                    add_3cycle_row_with_ub(ij, jk, ik);
                } else {
                    add_3cycle_row_with_lb(ij, jk, ik);
                }

                if (nof_new_rows >= PACE2024_CONST_NOF_CYCLE_CONSTRAINTS) {
                    return nof_new_rows;
                }
            }
        }

        return nof_new_rows;
    }

    //
    // 3-cycle helper methods
    //

    inline void check_3cycle(const int &i, const int &j, const int &k) {
        const int ij = get_variable_index(i, j);
        const int jk = get_variable_index(j, k);
        const int ik = get_variable_index(i, k);
        assert(ij < ik);
        assert(ik < jk);

        const double x = get_3cycle_value(ij, jk, ik);
        if (is_3cycle_lb_violated(x)) {
            get_bucket(-x).emplace_back(ij, jk, ik, false);
        }
        if (is_3cycle_ub_violated(x)) {
            get_bucket(x - 1).emplace_back(ij, jk, ik, true);
        }
    }

    bool check_3cycles() {
        PACE2024_DEBUG_PRINTF("\tstart check_3cycles\n");

        clear_buckets();
        bool break_for_loops;
        for (int i = 0; i < n1; ++i) {
            for (int j = i + 1; j < n1; ++j) {
                for (int k = j + 1; k < n1; ++k) {
                    assert(i < j);
                    assert(j < k);
                    check_3cycle(i, j, k);
                }
                break_for_loops = is_last_bucket_full();
                if (break_for_loops) break;
            }
            if (break_for_loops) break;
        }

        const int nof_new_rows = add_3cycle_rows();
        PACE2024_DEBUG_PRINTF("\tend   check_3cycles, number of new rows=%lld\n", nof_new_rows);
        return nof_new_rows > 0;
    }

    /**
     * @brief returns x_ij + x_jk - x_ik of the lp
     *
     * @param ij column index
     * @param jk column index
     * @param ik column index
     * @return double x_ij + x_jk - x_ik
     */
    inline double get_3cycle_value(const int &ij, const int &jk, const int &ik) {
        const double x_ij = glp_get_col_prim(lp, ij);
        const double x_jk = glp_get_col_prim(lp, jk);
        const double x_ik = glp_get_col_prim(lp, ik);
        return x_ij + x_jk - x_ik;
    }

    /**
     * @brief returns value of x < -1e-7
     */
    inline bool is_3cycle_lb_violated(const double &x) {
        return x < -params.tol_bnd;
    }

    /**
     * @brief returns value of x > 1 + 1e-7
     */
    inline bool is_3cycle_ub_violated(const double &x) {
        return x > 1. + params.tol_bnd;
    }

    //
    // row removal methods
    //

    /**
     * @brief returns lb != -DBL_MAX
     *
     * @param i row index
     * @return lb != -DBL_MAX
     */
    inline bool has_row_lb(const int &i) {
        const double lb = glp_get_row_lb(lp, i);
        return lb != -DBL_MAX;
    }

    /**
     * @brief returns ub != DBL_MAX
     *
     * @param i row index
     * @return ub != DBL_MAX
     */
    inline bool has_row_ub(const int &i) {
        const double ub = glp_get_row_ub(lp, i);
        return ub != DBL_MAX;
    }

    /**
     * @brief returns has_row_lb(i) && x > 0.
     *
     * @param i row index
     * @return has_row_lb(i) && x > 0.
     */
    inline bool has_row_lb_slack(const int &i) {
        const double x = glp_get_row_prim(lp, i);
        return has_row_lb(i) && x > 0.;
    }

    /**
     * @brief returns has_row_ub(i) && x < 1.
     *
     * @param i row index
     * @return has_row_ub(i) && x < 1.
     */
    inline bool has_row_ub_slack(const int &i) {
        const double x = glp_get_row_prim(lp, i);
        return has_row_ub(i) && x < 1.;
    }

    /**
     * @brief delete rows with positive slack from the lp
     * (each row just has either a lower or an upper bound)
     * (this is only slightly inefficient, since we delete positive slack rows, but it makes life easier)
     */
    inline void delete_positive_slack_rows() {
        PACE2024_DEBUG_PRINTF("\tstart delete_positive_slack_rows\n");

        std::vector<int> rows_to_remove;
        rows_to_remove.emplace_back(0);  // padding

        // gather rows to delete
        const int nof_rows = glp_get_num_rows(lp);
        for (int i = 1; i <= nof_rows; ++i) {
            if (has_row_lb_slack(i)) {
                rows_to_remove.emplace_back(i);
            }
            if (has_row_ub_slack(i)) {
                rows_to_remove.emplace_back(i);
            }
        }

        const int nof_removed_rows = static_cast<int>(rows_to_remove.size()) - 1;
        if (nof_removed_rows > 0) {
            PACE2024_DEBUG_PRINTF("\tend   delete_positive_slack_rows, number of removed rows=%lld\n", nof_removed_rows);
            glp_del_rows(lp, nof_removed_rows, &rows_to_remove[0]);
        }
    }

    //
    // integrality testing methods
    //

    /**
     * @brief checks if a column variable of the lp is integral
     *
     * @param j column index
     * @return true if integral (that is, it is in the params.tol_bnd open neighborhood of an integer)
     * @return false otherwise
     */
    inline bool
    is_column_integral(const int &j) {
        assert(0 <= j);
        assert(j < n1_choose_2);
        const double x = glp_get_col_prim(lp, j);
        constexpr double ub = 1. - PACE2024_CONST_INTEGER_TOLERANCE;
        if (x > PACE2024_CONST_INTEGER_TOLERANCE && x < ub) {
            return false;
        }
        return true;
    }

    //
    // helper methods
    //

    /**
     * @brief converts a vertex pair (i, j), i < j, to the lp column index
     *
     * @param i vertex
     * @param j vertex, i < j
     * @return int column index
     */
    inline int get_variable_index(const int &i, const int &j) {
        assert(i < j);
        int offset = n1_choose_2 - (n1 - i) * (n1 - i - 1) / 2;
        int index = offset + j - i - 1;
        assert(1 <= index && index <= n1_choose_2);
        return index;
    }
};

class highs_wrapper : public lp_wrapper {
    /// @brief lp solver
    Highs lp;
    /// @brief parameters for glp_simplex
    glp_smcp params;

   public:
    template <typename T>
    highs_wrapper(const general_bipartite_graph<T> &graph)
        : lp_wrapper(static_cast<int>(graph.get_n1())) {
        lp.setOptionValue("presolve", "off");
        lp.setOptionValue("solver", "simplex");
        lp.setOptionValue("parallel", "off");
        lp.setOptionValue("threads", 1);
        lp.setOptionValue("log_to_console", false);
        add_columns(graph);
    }

    highs_wrapper(const highs_wrapper &rhs) = delete;
    highs_wrapper(highs_wrapper &&rhs) = delete;
    highs_wrapper &operator=(const highs_wrapper &rhs) = delete;
    highs_wrapper &operator=(highs_wrapper &&rhs) = delete;

    ~highs_wrapper() {}

    /**
     * @brief search for violated constraints and add these to the lp
     *
     * @return true if successful
     * @return false otherwise
     */
    virtual bool cut() {
        bool success = check_3cycles();
        // todo: k-fence
        return success;
    }

    /**
     * @brief fixes column j to fix_to
     *
     * @param j column index
     * @param fix_to
     */
    virtual void fix_column(const int j, const double fix_to) {
        assert(0 <= j);
        assert(j < lp.getNumCol());
        PACE2024_DEBUG_PRINTF("fixed variable %5d to %1.0f\n", j, fix_to);
        lp.changeColBounds(j, fix_to, fix_to);
    }

    /**
     * @brief fixes columns based on the following:
     * if lower_bound = sum min(c_ij, c_ji), where c_ij and c_ji are crossing numbers,
     * and coeff = c_ij - c_ji is the objective coefficient of variable x_ij,
     * and abs(c_ij - c_ji) >= upper_bound - lower_bound,
     * then we are able to fix x_ij to 0 if coeff > 0 and to 1 otherwise
     *
     * @param new_upper_bound new upper bound on the optimal value
     */
    virtual void fix_columns(const int new_upper_bound) {
        assert(new_upper_bound < upper_bound);
        upper_bound = new_upper_bound;
        for (int j = 0; j < n1_choose_2; ++j) {
            const int diff = upper_bound - lower_bound;
            const double coeff = get_objective_coefficient(j);

            if (std::abs(coeff) >= static_cast<double>(diff) && coeff != 0.) {
                fix_column(j, coeff > 0 ? 0. : 1.);
                PACE2024_DEBUG_PRINTF("(optimality condition)\n");
            }
        }
    }

    /// @brief returns value of column j
    virtual double get_column_value(const int j) {
        if (j >= n1_choose_2) {
            PACE2024_DEBUG_PRINTF("now\n");
        }
        assert(0 <= j);
        assert(j < n1_choose_2);
        const HighsSolution &sol = lp.getSolution();
        return sol.col_value[j];
    }

    /// @brief returns the number of rows
    virtual int get_nof_cols() {
        return lp.getNumCol();
    }

    /// @brief returns the number of rows
    virtual int get_nof_rows() {
        return lp.getNumRow();
    }

    /// @brief returns the objective value of the lp
    virtual double get_objective_value() {
        return lp.getInfo().objective_function_value;
    }

    /// @brief returns the objective value of the lp rounded to the next integer
    virtual int64_t get_rounded_objective_value() {
        const double value = get_objective_value();
        assert(value >= 0);
        return llround(value);
    }

    /**
     * @brief returns a value based on the integrality of the current solution
     *
     * @return int 0, if solution is integral
     * @return int j, 1 <= j <= n1_choose_2, of nonintegral column otherwise
     */
    virtual int is_integral() {
        for (int j = 0; j < n1_choose_2; ++j) {
            if (!is_column_integral(j)) {
                return j;
            }
        }
        return 0;
    }

    /// @brief returns if an optimal feasible solution has been found
    virtual bool is_optimal() {
        const HighsModelStatus &model_status = lp.getModelStatus();
        return model_status == HighsModelStatus::kOptimal;
    }

    virtual void solve(bool delete_rows_after) {
        PACE2024_DEBUG_PRINTF("start glp_simplex\n");
        lp.run();
        PACE2024_DEBUG_PRINTF("end   glp_simplex, objective value=%f\n", get_objective_value());

        if (delete_rows_after) {
            delete_positive_slack_rows();
        }
    }

    /// @brief resets bounds of column j to 0 <= . <= 1
    virtual void unfix_column(const int j) {
        assert(0 <= j);
        assert(j < n1_choose_2);
        lp.changeColBounds(j, 0., 1.);
    }

   private:
    //
    // initialization methods
    //

    /**
     * @brief adds all columns of the ilp formulation of one-sided crossing minimization
     * and computes their objective coefficients (that is, the crossing numbers)
     *
     * @tparam T vertex param
     * @param graph the input graph
     */
    template <typename T>
    inline void add_columns(const general_bipartite_graph<T> &graph) {
        int obj_val_offset = 0;
        std::vector<double> lower_bounds(n1_choose_2, 0.);
        std::vector<double> upper_bounds(n1_choose_2, 1.);
        lp.addVars(n1_choose_2, &lower_bounds[0], &upper_bounds[0]);

        int k = 0;
        for (T i = 0; i < static_cast<T>(n1); ++i) {
            for (T j = i + 1; j < static_cast<T>(n1); ++j) {
                obj_val_offset += add_column(graph, i, j, k);
                ++k;
            }
        }
        // set constant term (shift/offset) in the objective function
        lp.changeObjectiveOffset(static_cast<double>(obj_val_offset));

        assert(n1_choose_2 == get_nof_cols());
    }

    /**
     * @brief computes the crossing numbers c_ij and c_ji,
     * sets the corresponding bounds for the column (perhaps even fixes them),
     *
     * @tparam T vertex type
     * @param graph the input graph
     * @param i vertex of free layer
     * @param j vertex of free layer
     * @param k column index
     * @return crossing number c_ji (for the objective offset)
     */
    template <typename T>
    inline int add_column(const general_bipartite_graph<T> &graph,
                          const T &i,
                          const T &j,
                          const int &k) {
        assert(k <= n1_choose_2);
        const auto [c_ij, c_ji] = crossing_numbers_of<T, int>(graph, i, j);

        // lower_bound = sum min(c_ij, c_ji)
        if (c_ij < c_ji) {
            lower_bound += c_ij;
        } else {
            lower_bound += c_ji;
        }

        if (c_ij == 0 && c_ji != 0) {
            // fix i < j in the ordering
            fix_column(k, 1.);
        } else if (c_ji == 0 && c_ij != 0) {
            // fix j < i in the ordering
            fix_column(k, 0.);
        } else {
            // set 0 <= x_ij <= 1
            unfix_column(k);
        }

        // todo: if d(i)=d(j), check if fixing is possible

        // set coefficient of added column
        if (c_ij != c_ji) {
            const double coeff = static_cast<double>(c_ij) - static_cast<double>(c_ji);
            lp.changeColCost(k, coeff);
        }

        return c_ji;
    }

    //
    // row addition methods
    //

    /**
     * @brief expects the violated 3-cycle ieqs in buckets.
     * adds the most violated <= PACE2024_CONST_NOF_CYCLE_CONSTRAINTS to the lp.
     *
     * @return int number of new rows
     */
    inline int add_3cycle_rows() {
        const int nof_new_rows = get_nof_bucket_entries();
        constexpr HighsInt array_length = PACE2024_CONST_NOF_CYCLE_CONSTRAINTS;
        constexpr HighsInt array_length_x3 = 3 * array_length;
        double lower_bounds[array_length];
        double upper_bounds[array_length];
        HighsInt starts[array_length];
        HighsInt indices[array_length_x3];
        double values[array_length_x3];

        int i = 0;
        bool go_on = true;
        for (auto r_it = buckets.rbegin(); r_it != buckets.rend() && go_on; ++r_it) {
            for (const auto &[ij, jk, ik, ub] : *r_it) {
                int i3 = 3 * i;
                (void)ub;  // suppress unused warning

                lower_bounds[i] = 0.;
                upper_bounds[i] = 1.;
                starts[i] = i3;
                indices[i3] = ij;
                indices[i3 + 1] = ik;
                indices[i3 + 2] = jk;
                values[i3] = 1.;
                values[i3 + 1] = -1.;
                values[i3 + 2] = 1.;

                ++i;
                if (i >= array_length) {
                    go_on = false;
                    break;
                }
            }
        }

        lp.addRows(nof_new_rows, lower_bounds, upper_bounds, 3 * nof_new_rows, starts, indices, values);
        return nof_new_rows;
    }

    //
    // 3-cycle helper methods
    //

    inline void check_3cycle(const int &i, const int &j, const int &k) {
        const int ij = get_variable_index(i, j);
        const int jk = get_variable_index(j, k);
        const int ik = get_variable_index(i, k);
        assert(ij < ik);
        assert(ik < jk);

        const double x = get_3cycle_value(ij, jk, ik);
        if (is_3cycle_lb_violated(x)) {
            get_bucket(-x).emplace_back(ij, jk, ik, false);
        }
        if (is_3cycle_ub_violated(x)) {
            get_bucket(x - 1).emplace_back(ij, jk, ik, true);
        }
    }

    bool check_3cycles() {
        PACE2024_DEBUG_PRINTF("\tstart check_3cycles\n");

        clear_buckets();
        bool break_for_loops;
        for (int i = 0; i < n1; ++i) {
            for (int j = i + 1; j < n1; ++j) {
                for (int k = j + 1; k < n1; ++k) {
                    assert(i < j);
                    assert(j < k);
                    check_3cycle(i, j, k);
                }
                break_for_loops = is_last_bucket_full();
                if (break_for_loops) break;
            }
            if (break_for_loops) break;
        }

        const int nof_new_rows = add_3cycle_rows();
        PACE2024_DEBUG_PRINTF("\tend   check_3cycles, number of new rows=%lld\n", nof_new_rows);
        return nof_new_rows > 0;
    }

    /**
     * @brief returns x_ij + x_jk - x_ik of the lp
     *
     * @param ij column index
     * @param jk column index
     * @param ik column index
     * @return double x_ij + x_jk - x_ik
     */
    inline double get_3cycle_value(const int &ij, const int &jk, const int &ik) {
        const double x_ij = get_column_value(ij);
        const double x_jk = get_column_value(jk);
        const double x_ik = get_column_value(ik);
        return x_ij + x_jk - x_ik;
    }

    /**
     * @brief returns value of x < -1e-7
     */
    inline bool is_3cycle_lb_violated(const double &x) {
        return x < -params.tol_bnd;
    }

    /**
     * @brief returns value of x > 1 + 1e-7
     */
    inline bool is_3cycle_ub_violated(const double &x) {
        return x > 1. + params.tol_bnd;
    }

    //
    // row removal methods
    //

    /**
     * @brief returns has_row_lb(i) && x > 0.
     *
     * @param i row index
     * @return has_row_lb(i) && x > 0.
     */
    inline bool has_row_slack(const int &i) {
        const HighsSolution &sol = lp.getSolution();
        const double x = sol.row_value[i];
        constexpr double ub = 1. - PACE2024_CONST_FEASIBILITY_TOLERANCE;
        return x > PACE2024_CONST_FEASIBILITY_TOLERANCE && x < ub;
    }

    /**
     * @brief delete rows with positive slack from the lp
     * (each row just has either a lower or an upper bound)
     * (this is only slightly inefficient, since we delete positive slack rows, but it makes life easier)
     */
    inline void delete_positive_slack_rows() {
        PACE2024_DEBUG_PRINTF("\tstart delete_positive_slack_rows\n");

        std::vector<int> rows_to_delete;

        // gather rows to delete
        const int nof_rows = get_nof_rows();
        for (int i = 0; i < nof_rows; ++i) {
            if (has_row_slack(i)) {
                rows_to_delete.emplace_back(i);
            }
        }

        const int nof_rows_to_delete = static_cast<int>(rows_to_delete.size());
        if (nof_rows_to_delete > 0) {
            PACE2024_DEBUG_PRINTF("\tend   delete_positive_slack_rows, number of removed rows=%d\n", nof_rows_to_delete);
            lp.deleteRows(nof_rows_to_delete, &rows_to_delete[0]);
        }
    }

    //
    // integrality testing methods
    //

    /**
     * @brief checks if a column variable of the lp is integral
     *
     * @param j column index
     * @return true if integral (that is, it is in the params.tol_bnd open neighborhood of an integer)
     * @return false otherwise
     */
    inline bool
    is_column_integral(const int &j) {
        const double x = get_column_value(j);
        constexpr double ub = 1. - PACE2024_CONST_INTEGER_TOLERANCE;
        if (x > PACE2024_CONST_INTEGER_TOLERANCE && x < ub) {
            return false;
        }
        return true;
    }

    //
    // helper methods
    //

    /**
     * @brief converts a vertex pair (i, j), i < j, to the lp column index
     *
     * @param i vertex
     * @param j vertex, i < j
     * @return int column index
     */
    inline int get_variable_index(const int &i, const int &j) {
        assert(i < j);
        int offset = n1_choose_2 - (n1 - i) * (n1 - i - 1) / 2;
        int index = offset + j - i - 1;
        assert(0 <= index && index < n1_choose_2);
        return index;
    }

    /// @brief returns the objective coefficient of column j
    inline double get_objective_coefficient(const int &j) {
        assert(0 <= j);
        assert(j < n1_choose_2);
        return lp.getLp().col_cost_[j];
    }
};

};  // namespace pace2024

#endif
