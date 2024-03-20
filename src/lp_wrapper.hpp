#ifndef PACE2024_LP_WRAPPER_HPP
#define PACE2024_LP_WRAPPER_HPP

#include <cassert>
#include <unordered_map>

#include "Highs.h"
#include "bipartite_graph.hpp"
#include "crossings.hpp"
#include "debug_printf.hpp"

#ifndef PACE2024_CONST_NOF_CYCLE_CONSTRAINTS
#define PACE2024_CONST_NOF_CYCLE_CONSTRAINTS 1024
#endif

#ifndef PACE2024_CONST_NOF_BUCKETS
#define PACE2024_CONST_NOF_BUCKETS 8
#endif

#ifndef PACE2024_CONST_FEASIBILITY_TOLERANCE
#define PACE2024_CONST_FEASIBILITY_TOLERANCE 1e-7
#endif

#ifndef PACE2024_CONST_INTEGER_TOLERANCE
#define PACE2024_CONST_INTEGER_TOLERANCE 1e-6
#endif

namespace pace2024 {

/**
 * @brief abstract wrapper class for an lp solver used by branch_and_cut
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
        assert(0 <= val);
        assert(val < 1);
        const int i = static_cast<int>(val * PACE2024_CONST_NOF_BUCKETS);
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

    // purely virtual methods
    virtual bool cut() = 0;
    virtual void delete_positive_slack_rows() = 0;
    virtual void fix_column(const int j, const double fix_to) = 0;
    virtual void fix_columns(const int new_upper_bound) = 0;
    virtual double get_variable_value(const int j) = 0;
    virtual int get_nof_cols() = 0;
    virtual int get_nof_rows() = 0;
    virtual double get_objective_value() = 0;
    virtual long get_rounded_objective_value() = 0;
    virtual inline bool is_column_integral(const int &j) = 0;
    virtual int is_integral() = 0;
    virtual bool is_optimal() = 0;
    virtual void solve() = 0;
    virtual void unfix_column(const int j) = 0;
};

class highs_wrapper : public lp_wrapper {
    /// @brief interface to lp model and solver
    Highs lp;

    /// @brief highs status field
    HighsStatus status;

    /// @brief number of rows before `cut()` added new rows
    int nof_old_rows{0};

    std::vector<int> magic;

    // for `delete_positive_slack_rows()`
    std::vector<HighsInt> rows_to_delete;

    // for `add_3cycle_rows()`
    std::vector<double> lower_bounds;
    std::vector<double> upper_bounds;
    std::vector<HighsInt> starts;
    std::vector<HighsInt> indices;
    std::vector<double> values;

   public:
    template <typename T>
    highs_wrapper(const bipartite_graph<T> &graph)
        : lp_wrapper(static_cast<int>(graph.get_n_free())),
          magic(n1_choose_2) {
        // configure highs lp solver
        status = lp.setOptionValue("presolve", "off");
        assert(status == HighsStatus::kOk);
        status = lp.setOptionValue("solver", "simplex");
        assert(status == HighsStatus::kOk);
        status = lp.setOptionValue("parallel", "off");
        assert(status == HighsStatus::kOk);
        status = lp.setOptionValue("threads", 1);
        assert(status == HighsStatus::kOk);
        status = lp.setOptionValue("log_to_console", false);
        assert(status == HighsStatus::kOk);

        add_variables(graph);

        rows_to_delete.reserve(PACE2024_CONST_NOF_CYCLE_CONSTRAINTS);
        lower_bounds.reserve(PACE2024_CONST_NOF_CYCLE_CONSTRAINTS);
        upper_bounds.reserve(PACE2024_CONST_NOF_CYCLE_CONSTRAINTS);
        starts.reserve(PACE2024_CONST_NOF_CYCLE_CONSTRAINTS);
        constexpr std::size_t PACE2024_CONST_NOF_CYCLE_CONSTRAINTS_x3 = 3 * PACE2024_CONST_NOF_CYCLE_CONSTRAINTS;
        indices.reserve(PACE2024_CONST_NOF_CYCLE_CONSTRAINTS_x3);
        values.reserve(PACE2024_CONST_NOF_CYCLE_CONSTRAINTS_x3);
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
        nof_old_rows = get_nof_rows();
        bool success = check_3cycles();
        // todo: k-fence
        return success;
    }

    /**
     * @brief delete rows with positive slack from the lp
     */
    virtual void delete_positive_slack_rows() {
        PACE2024_DEBUG_PRINTF("\tstart delete_positive_slack_rows\n");

        // gather rows to delete
        rows_to_delete.clear();
        for (int i = 0; i < nof_old_rows; ++i) {
            if (has_row_slack(i)) {
                rows_to_delete.emplace_back(i);
            }
        }

        const HighsInt nof_rows_to_delete = rows_to_delete.size();
        if (nof_rows_to_delete > 0) {
            lp.deleteRows(nof_rows_to_delete, &rows_to_delete[0]);
        }
        PACE2024_DEBUG_PRINTF("\tend   delete_positive_slack_rows, number of removed rows=%d\n", nof_rows_to_delete);
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
        assert(lower_bound <= new_upper_bound);
        upper_bound = new_upper_bound;
        for (int j = 0; j < lp.getNumCol(); ++j) {
            const int diff = upper_bound - lower_bound;
            const double coeff = get_objective_coefficient(j);

            if (std::abs(coeff) >= static_cast<double>(diff) && coeff != 0.) {
                fix_column(j, coeff > 0 ? 0. : 1.);
                PACE2024_DEBUG_PRINTF("(optimality condition)\n");
            }
        }
    }

    /// @brief returns value of column j
    virtual double get_variable_value(const int j) {
        assert(0 <= j);
        assert(j < n1_choose_2);
        if (magic[j] < 0) {
            return static_cast<double>(~magic[j]);
        } else {
            return lp.getSolution().col_value[magic[j]];
        }
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
    virtual long get_rounded_objective_value() {
        const double value = get_objective_value();
        assert(value >= 0);
        return lround(value);
    }

    /**
     * @brief returns a value based on the integrality of the current solution
     *
     * @return int -1, if solution is integral
     * @return int j, 0 <= j < lp.getNumCol(), of nonintegral column otherwise
     */
    virtual int is_integral() {
        for (int j = 0; j < lp.getNumCol(); ++j) {
            if (!is_column_integral(j)) {
                return j;
            }
        }
        return -1;
    }

    /// @brief returns if an optimal feasible solution has been found
    virtual bool is_optimal() {
        const HighsModelStatus &model_status = lp.getModelStatus();
        return model_status == HighsModelStatus::kOptimal;
    }

    /// @brief solves the lp
    virtual void solve() {
        PACE2024_DEBUG_PRINTF("start highs_simplex\n");
        lp.run();
        PACE2024_DEBUG_PRINTF("end   highs_simplex, objective value=%f, number of rows=%d\n",
                              get_objective_value(),
                              get_nof_rows());
    }

    /// @brief resets bounds of column j to 0 <= . <= 1
    virtual void unfix_column(const int j) {
        assert(0 <= j);
        assert(j < n1_choose_2);
        change_column_bounds(j, 0., 1.);
    }

   private:
    inline void change_column_bounds(const int &j, const double &&lb, const double &&ub) {
        status = lp.changeColBounds(j, lb, ub);
        assert(status == HighsStatus::kOk);
    }

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
    inline void add_variables(const bipartite_graph<T> &graph) {
        int obj_val_offset = 0;

        int k = 0;
        for (T i = 0; i < static_cast<T>(n1); ++i) {
            for (T j = i + 1; j < static_cast<T>(n1); ++j) {
                obj_val_offset += add_variable(graph, i, j, k);
                ++k;
            }
        }
        // set constant term (shift/offset) in the objective function
        lp.changeObjectiveOffset(static_cast<double>(obj_val_offset));

        assert(get_nof_cols() <= n1_choose_2);
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
    inline int add_variable(const bipartite_graph<T> &graph,
                            const T &i,
                            const T &j,
                            const int &k) {
        assert(k < n1_choose_2);
        const auto [c_ij, c_ji] = crossing_numbers_of(graph, i, j);

        // lower_bound = sum min(c_ij, c_ji)
        if (c_ij < c_ji) {
            lower_bound += c_ij;
        } else {
            lower_bound += c_ji;
        }

        // todo: if d(i)=d(j), check if fixing is possible
        if (c_ij == 0 && c_ji != 0) {
            // fix i < j in the ordering
            magic[k] = -2;
            return 0;
        } else if (c_ji == 0 && c_ij != 0) {
            // fix j < i in the ordering
            magic[k] = -1;
            return 0;
        } else {
            // set 0 <= x_ij <= 1
            const int l = lp.getNumCol();
            magic[k] = l;
            add_var();
            if (c_ij != c_ji) {
                change_column_cost(l, static_cast<double>(c_ij) - static_cast<double>(c_ji));
            }
            return c_ji;
        }
    }

    inline void add_var() {
        status = lp.addVar(0., 1.);
        assert(status == HighsStatus::kOk);
    }

    inline void change_column_cost(const int &l, const double &&cost) {
        status = lp.changeColCost(l, cost);
        assert(status == HighsStatus::kOk);
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
        const HighsInt nof_new_rows = std::min(get_nof_bucket_entries(), PACE2024_CONST_NOF_CYCLE_CONSTRAINTS);
        if (nof_new_rows <= 0) return 0;

        lower_bounds.clear();
        upper_bounds.clear();
        starts.clear();
        indices.clear();
        values.clear();

        int i = 0;
        bool go_on = true;
        for (auto r_it = buckets.rbegin(); r_it != buckets.rend() && go_on; ++r_it) {
            for (const auto &[ij, jk, ik, ub] : *r_it) {
                (void)ub;  // suppress unused warning

                const double bound_offset = get_3cycle_bound_offset(ij, jk, ik);
                lower_bounds.emplace_back(0. + bound_offset);
                upper_bounds.emplace_back(1. + bound_offset);
                starts.emplace_back(indices.size());
                if (magic[ij] >= 0) {
                    indices.emplace_back(magic[ij]);
                    values.emplace_back(1.);
                }
                if (magic[ik] >= 0) {
                    indices.emplace_back(magic[ik]);
                    values.emplace_back(-1.);
                }
                if (magic[jk] >= 0) {
                    indices.emplace_back(magic[jk]);
                    values.emplace_back(1.);
                }
                assert(magic[ij] >= 0 || magic[ik] >= 0 || magic[jk] >= 0);

                ++i;
                if (i >= PACE2024_CONST_NOF_CYCLE_CONSTRAINTS) {
                    go_on = false;
                    break;
                }
            }
        }
        assert(i == nof_new_rows);

        lp.addRows(nof_new_rows,
                   &lower_bounds[0], &upper_bounds[0],
                   indices.size(),
                   &starts[0], &indices[0], &values[0]);
        return nof_new_rows;
    }

    //
    // 3-cycle helper methods
    //

    /**
     * @brief checks if the 3-cycle inequality for ijk/ikj is violated
     *
     * @return true if so
     * @return false otherwise
     */
    inline bool check_3cycle(const int &i, const int &j, const int &k) {
        const int ij = get_variable_index(i, j);
        const int jk = get_variable_index(j, k);
        const int ik = get_variable_index(i, k);
        assert(ij < ik);
        assert(ik < jk);

        const double x = get_3cycle_value(ij, jk, ik);
        constexpr double interval_width = 1. + 2e-7;
        if (is_3cycle_lb_violated(x)) {
            const double x_normalized = -x / interval_width;
            get_bucket(x_normalized).emplace_back(ij, jk, ik, false);
            return true;
        }
        if (is_3cycle_ub_violated(x)) {
            constexpr double ub = 1. + PACE2024_CONST_FEASIBILITY_TOLERANCE;
            const double x_normalized = (x - ub) / interval_width;
            get_bucket(x_normalized).emplace_back(ij, jk, ik, true);
            return true;
        }
        return false;
    }

    /**
     * @brief checks if the 3-cycle inequalities are violated
     *
     * @return true if so
     * @return false otherwise
     */
    bool check_3cycles() {
        PACE2024_DEBUG_PRINTF("\tstart check_3cycles\n");

        clear_buckets();
        bool break_for_loops;
        for (int i = 0; i < n1 - 2; ++i) {
            for (int j = i + 1; j < n1 - 1; ++j) {
                for (int k = j + 1; k < n1; ++k) {
                    assert(i < j);
                    assert(j < k);
                    bool violated = check_3cycle(i, j, k);
                    if (violated) {
                        break_for_loops = is_last_bucket_full();
                        if (break_for_loops) {
                            PACE2024_DEBUG_PRINTF("\t\tlast bucket full\n");
                            break;
                        }
                    }
                }
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
        const double x_ij = get_variable_value(ij);
        const double x_jk = get_variable_value(jk);
        const double x_ik = get_variable_value(ik);
        return x_ij + x_jk - x_ik;
    }

    /**
     * @brief returns value of x < -1e-7
     */
    inline bool is_3cycle_lb_violated(const double &x) {
        return x < -PACE2024_CONST_FEASIBILITY_TOLERANCE;
    }

    /**
     * @brief returns value of x > 1 + 1e-7
     */
    inline bool is_3cycle_ub_violated(const double &x) {
        constexpr double ub = 1. + PACE2024_CONST_FEASIBILITY_TOLERANCE;
        return x > ub;
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
        const double x = get_row_value(i);
        const auto [lb, ub] = get_row_bounds(i);
        return x > lb + PACE2024_CONST_FEASIBILITY_TOLERANCE && x < ub - PACE2024_CONST_FEASIBILITY_TOLERANCE;
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
    inline bool is_column_integral(const int &j) {
        assert(0 <= j);
        assert(j < lp.getNumCol());
        const double x = lp.getSolution().col_value[j];
        constexpr double ub = 1. - PACE2024_CONST_INTEGER_TOLERANCE;
        if (x > PACE2024_CONST_INTEGER_TOLERANCE && x < ub) {
            return false;
        }
        return true;
    }

    //
    // helper methods
    //

    /// @brief returns the objective coefficient of column j
    inline double get_objective_coefficient(const int &j) {
        assert(0 <= j);
        assert(j < lp.getNumCol());
        return lp.getLp().col_cost_[j];
    }

    /**
     * @brief returns the lower and upper bound offset for 3-cycle inequalities ijk
     * due to permanently fixed variables that are not incorporated in the lp
     */
    inline int get_3cycle_bound_offset(const int &ij, const int &jk, const int &ik) {
        double offset = 0.;
        if (magic[ij] < 0) {
            offset -= ~magic[ij];
        }
        if (magic[ik] < 0) {
            offset += ~magic[ik];
        }
        if (magic[jk] < 0) {
            offset -= ~magic[jk];
        }
        return offset;
    }

    /**
     * @brief returns the lower and upper bound offset for 3-cycle inequality with index i
     * due to permanently fixed variables that are not incorporated in the lp
     */
    inline std::pair<double, double> get_row_bounds(const int &i) {
        assert(0 <= i);
        assert(i < lp.getNumRow());
        const HighsLp &lp_ = lp.getLp();
        const double lb = lp_.row_lower_[i];
        const double ub = lp_.row_upper_[i];
        return std::make_pair(lb, ub);
    }

    /// @brief returns the value of row i
    inline double get_row_value(const int &i) {
        return lp.getSolution().row_value[i];
    }
};

};  // namespace pace2024

#endif
