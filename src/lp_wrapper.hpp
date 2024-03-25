#ifndef PACE_LP_WRAPPER_HPP
#define PACE_LP_WRAPPER_HPP

#include <cassert>
#include <limits>

#include "Highs.h"
#include "debug_printf.hpp"
#include "index.hpp"
#include "instance.hpp"

#ifndef PACE_CONST_NOF_CYCLE_CONSTRAINTS
#define PACE_CONST_NOF_CYCLE_CONSTRAINTS 1024u
#endif

#ifndef PACE_CONST_NOF_BUCKETS
#define PACE_CONST_NOF_BUCKETS 8u
#endif

#ifndef PACE_CONST_FEASIBILITY_TOLERANCE
#define PACE_CONST_FEASIBILITY_TOLERANCE 1e-7
#endif

#ifndef PACE_CONST_INTEGER_TOLERANCE
#define PACE_CONST_INTEGER_TOLERANCE 1e-6
#endif

namespace pace {

/**
 * @brief abstract wrapper class for an lp solver used by branch_and_cut
 */
class lp_wrapper {
   protected:
    //
    // bucket attributes and methods
    //

    using bucket_entry = std::tuple<HighsInt, HighsInt, HighsInt>;

    /// @brief buckets for bucket sorting violated 3-cycle inequalities
    std::vector<std::vector<bucket_entry>> buckets{PACE_CONST_NOF_BUCKETS};

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
        const std::size_t i = static_cast<std::size_t>(val * PACE_CONST_NOF_BUCKETS);
        assert(i < PACE_CONST_NOF_BUCKETS);
        return buckets[i];
    }

    /// @brief returns sum of the number of elements in each bucket of `buckets`
    inline std::size_t get_nof_bucket_entries() {
        std::size_t n = 0;
        for (const auto &bucket : buckets) {
            n += bucket.size();
        }
        return n;
    }

    /**
     * @brief returns if the last bucket is full
     *
     * @return true if last bucket has >= PACE_CONST_NOF_CYCLE_CONSTRAINTS elements
     * @return false otherwise
     */
    inline bool is_last_bucket_full() {
        return (*buckets.rbegin()).size() >= PACE_CONST_NOF_CYCLE_CONSTRAINTS;
    }

    //
    // lp related attributes
    //

    //
    // instance related attributes
    //

    /// @brief number of vertices in free layer
    const std::size_t n1;

    /// @brief n1_choose_2 = n1 * (n1 - 1) / 2
    const std::size_t n1_choose_2;

    /// @brief lower bound
    const uint32_t &lower_bound;

    /// @brief best upper bound
    uint32_t upper_bound{std::numeric_limits<uint32_t>::max()};

   public:
    lp_wrapper(const std::size_t n1, const uint32_t &lower_bound)
        : n1(n1), n1_choose_2(n1 * (n1 - 1) / 2), lower_bound(lower_bound) {
        assert(n1 >= 2);
    }

    virtual ~lp_wrapper() {}

    // purely virtual methods
    virtual bool cut() = 0;
    virtual void delete_positive_slack_rows() = 0;
    virtual void fix_column(const HighsInt &j, const double fix_to) = 0;
    virtual void fix_columns(const uint32_t new_upper_bound) = 0;
    virtual double get_column_value(const HighsInt j) = 0;
    virtual std::size_t get_nof_cols() = 0;
    virtual std::size_t get_nof_rows() = 0;
    virtual double get_objective_value() = 0;
    virtual uint32_t get_rounded_objective_value() = 0;
    virtual double get_variable_value(const std::size_t &uv) = 0;
    virtual inline bool is_column_integral(const HighsInt &j) = 0;
    virtual HighsInt is_integral() = 0;
    virtual bool is_optimal() = 0;
    virtual void run() = 0;
    virtual void unfix_column(const HighsInt &j) = 0;
};

/**
 * @brief
 *
 * @tparam T vertex type
 */
template <typename T>
class highs_wrapper : public lp_wrapper {
    /// @brief interface to lp model and solver
    Highs lp;

    /// @brief highs status field
    HighsStatus status;

    /// @brief number of rows before `cut()` added new rows
    std::size_t nof_old_rows{0};

    /// @brief see branch_and_cut class
    const std::vector<int> &magic;

    // for `delete_positive_slack_rows()`
    std::vector<HighsInt> rows_to_delete;

    // for `add_3cycle_rows()`
    std::vector<double> lower_bounds;
    std::vector<double> upper_bounds;
    std::vector<HighsInt> starts;
    std::vector<HighsInt> indices;
    std::vector<double> values;

   public:
    template <typename R>
    highs_wrapper(const instance<T, R> &instance, const std::vector<int> &magic)
        : lp_wrapper(instance.graph().get_n_free(), instance.get_lower_bound()), magic(magic) {
        // configure highs lp solver
        status = lp.setOptionValue("presolve", "off");
        assert(status == HighsStatus::kOk);
        status = lp.setOptionValue("solver", "simplex");
        assert(status == HighsStatus::kOk);
        status = lp.setOptionValue("parallel", "off");
        assert(status == HighsStatus::kOk);
        status = lp.setOptionValue("log_to_console", false);
        assert(status == HighsStatus::kOk);

        rows_to_delete.reserve(PACE_CONST_NOF_CYCLE_CONSTRAINTS);
        lower_bounds.reserve(PACE_CONST_NOF_CYCLE_CONSTRAINTS);
        upper_bounds.reserve(PACE_CONST_NOF_CYCLE_CONSTRAINTS);
        starts.reserve(PACE_CONST_NOF_CYCLE_CONSTRAINTS);
        indices.reserve(3 * PACE_CONST_NOF_CYCLE_CONSTRAINTS);
        values.reserve(3 * PACE_CONST_NOF_CYCLE_CONSTRAINTS);
    }

    highs_wrapper(const highs_wrapper &rhs) = delete;
    highs_wrapper(highs_wrapper &&rhs) = delete;
    highs_wrapper &operator=(const highs_wrapper &rhs) = delete;
    highs_wrapper &operator=(highs_wrapper &&rhs) = delete;

    ~highs_wrapper() {}

    //
    // virtual methods
    //

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
        PACE_DEBUG_PRINTF("\tstart delete_positive_slack_rows\n");

        // gather rows to delete
        rows_to_delete.clear();
        for (std::size_t i = 0; i < nof_old_rows; ++i) {
            if (has_row_slack(i)) {
                rows_to_delete.emplace_back(i);
            }
        }

        const HighsInt nof_rows_to_delete = rows_to_delete.size();
        if (nof_rows_to_delete > 0) {
            lp.deleteRows(nof_rows_to_delete, &rows_to_delete[0]);
        }
        PACE_DEBUG_PRINTF("\tend   delete_positive_slack_rows, number of removed rows=%d\n",
                          nof_rows_to_delete);
    }

    /**
     * @brief fixes column j to fix_to
     *
     * @param j column index
     * @param fix_to value to assign to column j
     */
    virtual void fix_column(const HighsInt &j, const double fix_to) {
        assert(0 <= j);
        assert(j < lp.getNumCol());
        assert(0. <= fix_to);
        assert(fix_to <= 1.);
        PACE_DEBUG_PRINTF("fixed variable %5d to %1.0f\n", j, fix_to);
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
    virtual void fix_columns(const uint32_t new_upper_bound) {
        assert(new_upper_bound < upper_bound);
        assert(lower_bound <= new_upper_bound);
        upper_bound = new_upper_bound;
        for (HighsInt j = 0; j < lp.getNumCol(); ++j) {
            const uint32_t diff = upper_bound - lower_bound;
            const double coeff = get_objective_coefficient(j);

            if (std::abs(coeff) >= static_cast<double>(diff) && coeff != 0.) {
                fix_column(j, coeff > 0 ? 0. : 1.);
                PACE_DEBUG_PRINTF("(optimality condition)\n");
            }
        }
    }

    virtual double get_column_value(const HighsInt j) {
        assert(0 <= j);
        assert(j < lp.getNumCol());
        return lp.getSolution().col_value[j];
    }

    /// @brief returns value of variable x_uv
    virtual double get_variable_value(const std::size_t &uv) {
        assert(uv < n1_choose_2);
        if (magic[uv] < 0) {
            return static_cast<double>(~magic[uv]);
        } else {
            return lp.getSolution().col_value[magic[uv]];
        }
    }

    /// @brief returns the number of rows
    virtual std::size_t get_nof_cols() { return lp.getNumCol(); }

    /// @brief returns the number of rows
    virtual std::size_t get_nof_rows() { return lp.getNumRow(); }

    /// @brief returns the objective value of the lp
    virtual double get_objective_value() { return lp.getInfo().objective_function_value; }

    /// @brief returns the objective value of the lp rounded to the next integer
    virtual uint32_t get_rounded_objective_value() {
        const double value = get_objective_value();
        assert(value >= 0);
        return static_cast<uint32_t>(lround(value));
    }

    /**
     * @brief returns a value based on the integrality of the current solution
     *
     * @return int -1, if solution is integral
     * @return int j, 0 <= j < lp.getNumCol(), of nonintegral column otherwise
     */
    virtual HighsInt is_integral() {
        for (HighsInt j = 0; j < lp.getNumCol(); ++j) {
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
    virtual void run() {
        PACE_DEBUG_PRINTF("start highs_simplex\n");
        lp.run();
        PACE_DEBUG_PRINTF("end   highs_simplex\n");
    }

    /// @brief resets bounds of column j to 0 <= . <= 1
    virtual void unfix_column(const HighsInt &j) {
        assert(0 <= j);
        assert(static_cast<std::size_t>(j) < n1_choose_2);
        change_column_bounds(j, 0., 1.);
    }

   private:
    inline void change_column_bounds(const HighsInt &j, const double &&lb, const double &&ub) {
        status = lp.changeColBounds(j, lb, ub);
        assert(status == HighsStatus::kOk);
    }

    //
    // initialization methods
    //

    /**
     * @brief adds all variables of the ilp formulation of one-sided crossing minimization
     * and computes their objective coefficients from the crossing numbers
     *
     * @tparam R crossing number type
     * @param instance the instance
     */
    template <typename R>
    inline void add_columns(const instance<T, R> &instance) {
        double objetive_offset = 0.;
        const folded_matrix<R> &cr_matrix = instance.cr_matrix();

        // add variables which are not yet settled
        std::size_t k = 0;
        for (T u = 0; u < n1; ++u) {
            for (T v = u + 1; v < n1; ++v) {
                assert(k < n1_choose_2);

                if (magic[k] >= 0) {
                    const HighsInt l = lp.getNumCol();
                    assert(l == magic[k]);
                    add_column();

                    const double &c_uv = cr_matrix(u, v);
                    const double &c_vu = cr_matrix(v, u);
                    objetive_offset += c_vu;
                    if (c_uv != c_vu) {
                        change_column_cost(l, c_uv - c_vu);
                    }
                }

                ++k;
            }
        }

        // set constant term (shift/offset) in the objective function
        lp.changeObjectiveOffset(objetive_offset);

        assert(get_nof_cols() <= n1_choose_2);
    }

    inline void add_column() {
        status = lp.addVar(0., 1.);
        assert(status == HighsStatus::kOk);
    }

    inline void change_column_cost(const HighsInt &l, const double &&cost) {
        status = lp.changeColCost(l, cost);
        assert(status == HighsStatus::kOk);
    }

    //
    // row addition methods
    //

    /**
     * @brief expects the violated 3-cycle ieqs in buckets.
     * adds the most violated <= PACE_CONST_NOF_CYCLE_CONSTRAINTS to the lp.
     *
     * @return std::size_t number of new rows
     */
    inline std::size_t add_3cycle_rows() {
        const std::size_t nof_new_rows = std::min(
            get_nof_bucket_entries(), static_cast<std::size_t>(PACE_CONST_NOF_CYCLE_CONSTRAINTS));
        if (nof_new_rows <= 0) return 0;

        lower_bounds.clear();
        upper_bounds.clear();
        starts.clear();
        indices.clear();
        values.clear();

        std::size_t i = 0;
        bool go_on = true;
        for (auto r_it = buckets.rbegin(); r_it != buckets.rend() && go_on; ++r_it) {
            for (const auto &[ij, jk, ik] : *r_it) {
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
                if (i >= PACE_CONST_NOF_CYCLE_CONSTRAINTS) {
                    go_on = false;
                    break;
                }
            }
        }
        assert(i == nof_new_rows);

        lp.addRows(nof_new_rows, &lower_bounds[0], &upper_bounds[0],  //
                   indices.size(), &starts[0], &indices[0], &values[0]);
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
    inline bool check_3cycle(const T &u, const T &v, const T &w) {
        const std::size_t uv = flat_index(n1, n1_choose_2, u, v);
        const std::size_t vw = flat_index(n1, n1_choose_2, v, w);
        const std::size_t uw = flat_index(n1, n1_choose_2, u, w);
        assert(uv < uw);
        assert(uw < vw);

        const double x = get_3cycle_value(uv, vw, uw);
        constexpr double interval_width = 1. + 2e-7;
        if (is_3cycle_lb_violated(x)) {
            const double x_normalized = -x / interval_width;
            get_bucket(x_normalized).emplace_back(uv, vw, uw);
            return true;
        }
        if (is_3cycle_ub_violated(x)) {
            constexpr double ub = 1. + PACE_CONST_FEASIBILITY_TOLERANCE;
            const double x_normalized = (x - ub) / interval_width;
            get_bucket(x_normalized).emplace_back(uv, vw, uw);
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
        PACE_DEBUG_PRINTF("\tstart check_3cycles\n");

        clear_buckets();
        bool break_for_loops;
        for (T u = 0; u < n1 - 2; ++u) {
            for (T v = u + 1; v < n1 - 1; ++v) {
                for (T w = v + 1; w < n1; ++w) {
                    assert(u < v);
                    assert(v < w);
                    bool violated = check_3cycle(u, v, w);
                    if (violated) {
                        break_for_loops = is_last_bucket_full();
                        if (break_for_loops) {
                            PACE_DEBUG_PRINTF("\t\tlast bucket full\n");
                            break;
                        }
                    }
                }
                if (break_for_loops) break;
            }
            if (break_for_loops) break;
        }

        const int nof_new_rows = add_3cycle_rows();
        PACE_DEBUG_PRINTF("\tend   check_3cycles, number of new rows=%lld\n", nof_new_rows);
        return nof_new_rows > 0;
    }

    /**
     *  @brief returns x_uv + x_vw - x_uw of the lp
     */
    inline double get_3cycle_value(const int &uv, const int &vw, const int &uw) {
        const double x_uv = get_variable_value(uv);
        const double x_vw = get_variable_value(vw);
        const double x_uw = get_variable_value(uw);
        return x_uv + x_vw - x_uw;
    }

    /**
     * @brief returns value of x < -1e-7
     */
    inline bool is_3cycle_lb_violated(const double &x) {
        return x < -PACE_CONST_FEASIBILITY_TOLERANCE;
    }

    /**
     * @brief returns value of x > 1 + 1e-7
     */
    inline bool is_3cycle_ub_violated(const double &x) {
        constexpr double ub = 1. + PACE_CONST_FEASIBILITY_TOLERANCE;
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
    inline bool has_row_slack(const std::size_t &i) {
        const double x = get_row_value(i);
        const auto [lb, ub] = get_row_bounds(i);
        return x > lb + PACE_CONST_FEASIBILITY_TOLERANCE &&
               x < ub - PACE_CONST_FEASIBILITY_TOLERANCE;
    }

    //
    // integrality testing methods
    //

    /**
     * @brief checks if a column variable of the lp is integral
     *
     * @param j column index
     * @return true if integral (that is, it is in the params.tol_bnd open neighborhood of an
     * integer)
     * @return false otherwise
     */
    inline bool is_column_integral(const HighsInt &j) {
        assert(0 <= j);
        assert(j < lp.getNumCol());
        const double x = lp.getSolution().col_value[j];
        constexpr double ub = 1. - PACE_CONST_INTEGER_TOLERANCE;
        if (x > PACE_CONST_INTEGER_TOLERANCE && x < ub) {
            return false;
        }
        return true;
    }

    //
    // helper methods
    //

    /// @brief returns the objective coefficient of column j
    inline double get_objective_coefficient(const HighsInt &j) {
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
    inline double get_row_value(const int &i) { return lp.getSolution().row_value[i]; }
};

};  // namespace pace

#endif
