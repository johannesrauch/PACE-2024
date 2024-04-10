#ifndef PACE_LP_WRAPPER_HPP
#define PACE_LP_WRAPPER_HPP

#include <cassert>
#include <limits>
#include <unordered_map>

#include "Highs.h"
#include "debug_printf.hpp"
#include "index.hpp"
#include "instance.hpp"

#ifndef PACE_CONST_NOF_CYCLE_CONSTRAINTS
#define PACE_CONST_NOF_CYCLE_CONSTRAINTS 256u
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

template <typename T>
struct highs_wrapper_info {
    const T &u_old;
    const T &v_old;
    const T &w_old;

    std::size_t nof_rows{0};
    std::size_t nof_deleted_rows{0};
    std::size_t nof_delete_rows_spared{0};
    std::size_t nof_added_rows{0};
    std::size_t nof_iterations_simplex{0};
    std::size_t nof_iterations_3cycles{0};
    std::size_t nof_bucket_entries{0};

    double t_simplex{0.};
    double objective_value{0.};
};

/**
 * @brief
 *
 * @tparam T vertex type
 */
template <typename T>
class highs_wrapper {
    /// @brief number of vertices in free layer
    const std::size_t n;

    /// @brief n_choose_2 = n * (n - 1) / 2
    const std::size_t n_choose_2;

    /// @brief lower bound
    const uint32_t lower_bound;

    /// @brief best upper bound
    uint32_t upper_bound;

    /// @brief interface to lp model and solver
    Highs lp;

    /// @brief highs status field
    HighsStatus status{HighsStatus::kOk};

    /// @brief number of rows before cut() added new rows
    std::size_t nof_old_rows{0};

    /// @brief see branch_and_cut class
    const std::vector<int> &magic;

    // for delete_positive_slack_rows()
    std::vector<HighsInt> rows_to_delete;

    // for add_3cycle_rows()
    std::vector<double> lower_bounds;
    std::vector<double> upper_bounds;
    std::vector<HighsInt> starts;
    std::vector<HighsInt> indices;
    std::vector<double> values;

    T u_old{0}, v_old{1}, w_old{2};

    //
    // bucket attributes and methods
    //

    using bucket_entry = std::tuple<T, T, T>;

    /// @brief buckets for bucket sorting violated 3-cycle inequalities
    std::vector<std::vector<bucket_entry>> buckets{PACE_CONST_NOF_BUCKETS};

    //
    // info
    //

    highs_wrapper_info<T> info;

    //
    // triples
    //

    struct ordered_triple {
        T u, v, w;

        ordered_triple(const T &u, const T &v, const T &w) : u(u), v(v), w(w) {}

        bool operator==(const ordered_triple &other) const { return u == other.u && v == other.v && w == other.w; }
    };

    struct row_info {
        uint16_t nof_times_deleted{0};
        uint16_t nof_times_spared{0};
    };

    struct ordered_triple_hash {
        std::size_t operator()(const ordered_triple &t) const {
            const std::hash<T> hash;
            return hash(t.u) ^ hash(t.v) ^ hash(t.w);
        }
    };

    std::unordered_map<ordered_triple, row_info, ordered_triple_hash> rows_info;

    std::list<ordered_triple> rows;

    std::size_t nof_max_new_rows{PACE_CONST_NOF_CYCLE_CONSTRAINTS};

   public:
    template <typename R>
    highs_wrapper(const instance<T, R> &instance,                 //
                  const std::vector<int> &magic,                  //
                  const std::vector<std::pair<T, T>> &unsettled,  //
                  const uint32_t objective_offset,                //
                  const uint32_t upper_bound)
        : n(instance.graph().get_n_free()),
          n_choose_2(n * (n - 1) / 2),
          lower_bound(instance.get_lower_bound()),
          upper_bound(upper_bound),
          magic(magic),
          info{u_old, v_old, w_old} {
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

        add_columns(instance, unsettled, objective_offset);
    }

    highs_wrapper(const highs_wrapper &rhs) = delete;
    highs_wrapper(highs_wrapper &&rhs) = delete;
    highs_wrapper &operator=(const highs_wrapper &rhs) = delete;
    highs_wrapper &operator=(highs_wrapper &&rhs) = delete;

    ~highs_wrapper() {}

    //
    // lp solving methods
    //

    /**
     * @brief search for violated constraints and add these to the lp
     *
     * @return true if successful
     * @return false otherwise
     */
    bool cut() {
        nof_old_rows = get_nof_rows();
        bool success = check_3cycles();
        return success;
    }

    /// @brief solves the lp
    void run() {
        lp.run();
        info.nof_rows = get_nof_rows();
        info.nof_iterations_simplex = lp.getSimplexIterationCount();
        info.t_simplex = lp.getRunTime();
        info.objective_value = get_objective_value();
    }

    //
    // row modification methods
    //

    /**
     * @brief delete rows with positive slack from the lp
     */
    void delete_positive_slack_rows() {
        info.nof_delete_rows_spared = 0;

        // gather rows to delete
        rows_to_delete.clear();
        auto it_rows = rows.begin();
        for (std::size_t i = 0; i < nof_old_rows; ++i) {
            assert(it_rows != rows.end());

            const bool has_slack = has_row_slack(i);
            bool spare = false;
            if (has_slack) {
                auto it = rows_info.find(*it_rows);
                if (it == rows_info.end()) {
                    bool success;
                    std::tie(it, success) = rows_info.emplace(*it_rows, row_info{});
                    assert(success);
                }
                assert(it != rows_info.end());
                spare = it->second.nof_times_deleted > it->second.nof_times_spared;

                if (!spare) {
                    rows_to_delete.emplace_back(i);
                    ++it->second.nof_times_deleted;
                    it->second.nof_times_spared = 0;
                    it_rows = rows.erase(it_rows);
                } else {
                    ++it_rows;
                    ++it->second.nof_times_spared;
                    ++info.nof_delete_rows_spared;
                }
            } else {
                ++it_rows;
            }
        }

        const HighsInt nof_rows_to_delete = rows_to_delete.size();
        info.nof_deleted_rows = rows_to_delete.size();
        if (nof_rows_to_delete > 0) {
            lp.deleteRows(nof_rows_to_delete, &rows_to_delete[0]);
        }
    }

    //
    // column modification methods
    //

    /**
     * @brief adds all variables of the ilp formulation of one-sided crossing minimization
     * and computes their objective coefficients from the crossing numbers
     *
     * @tparam R crossing number type
     * @param instance the instance
     */
    template <typename R>
    inline void add_columns(const instance<T, R> &instance,                 //
                            const std::vector<std::pair<T, T>> &unsettled,  //
                            uint32_t objective_offset) {
        const folded_matrix<R> &cr_matrix = instance.cr_matrix();

        // add variables which are not yet settled
        for (const auto &[u, v] : unsettled) {
            const HighsInt l = lp.getNumCol();
            assert(l == magic[flat_index(n, n_choose_2, u, v)]);
            add_column();

            const R &c_uv = cr_matrix(u, v);
            const R &c_vu = cr_matrix(v, u);
            objective_offset += c_vu;
            if (c_uv != c_vu) {
                change_column_cost(l, static_cast<double>(c_uv) - static_cast<double>(c_vu));
            }
        }

        // set constant term (shift/offset) in the objective function
        lp.changeObjectiveOffset(objective_offset);

        assert(get_nof_cols() <= n_choose_2);
    }

    /// @brief adds a column with bounds 0. <= (*) <= 1.
    inline void add_column() {
        status = lp.addVar(0., 1.);
        assert(status == HighsStatus::kOk);
    }

    inline void change_column_bounds(const std::size_t j, const double lb, const double ub) {
        assert(j < get_nof_cols());
        status = lp.changeColBounds(j, lb, ub);
        assert(status == HighsStatus::kOk);
    }

    inline void change_column_cost(const std::size_t j, const double cost) {
        status = lp.changeColCost(j, cost);
        assert(status == HighsStatus::kOk);
    }

    /**
     * @brief fixes column j to fix_to
     *
     * @param j column index
     * @param fix_to value to assign to column j
     */
    void fix_column(const std::size_t j, const double fix_to) {
        assert(0. <= fix_to);
        assert(fix_to <= 1.);
        PACE_DEBUG_PRINTF("\tfixed variable %5d to %1.0f\n", j, fix_to);
        change_column_bounds(j, fix_to, fix_to);
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
    void fix_columns(const uint32_t new_upper_bound) {
        assert(new_upper_bound < upper_bound);
        assert(lower_bound <= new_upper_bound);
        upper_bound = new_upper_bound;
        for (std::size_t j = 0; j < get_nof_cols(); ++j) {
            const uint32_t diff = upper_bound - lower_bound;
            const double coeff = get_objective_coefficient(j);

            if (std::abs(coeff) >= static_cast<double>(diff) && coeff != 0.) {
                fix_column(j, coeff > 0 ? 0. : 1.);
                PACE_DEBUG_PRINTF("(optimality condition)\n");
            }
        }
    }

    /// @brief resets bounds of column j to 0 <= . <= 1
    void unfix_column(const std::size_t j) { change_column_bounds(j, 0., 1.); }

    //
    // getter
    //

    const highs_wrapper_info<T> &get_info() { return info; }

    void get_columns(std::vector<double> &col_value) const { col_value = lp.getSolution().col_value; }

    /// @brief returns value of column j
    double get_column_value(const std::size_t j) {
        assert(j < get_nof_cols());
        return lp.getSolution().col_value[j];
    }

    /// @brief returns value of variable x_uv
    double get_variable_value(const std::size_t uv) {
        assert(uv < n_choose_2);
        if (magic[uv] < 0) {
            return static_cast<double>(~magic[uv]);
        } else {
            return lp.getSolution().col_value[magic[uv]];
        }
    }

    /// @brief returns the number of rows
    std::size_t get_nof_cols() const { return lp.getNumCol(); }

    /// @brief returns the number of rows
    std::size_t get_nof_rows() const { return lp.getNumRow(); }

    /// @brief returns the objective value of the lp
    double get_objective_value() { return lp.getInfo().objective_function_value; }

    /// @brief returns the objective value of the lp rounded to the next integer
    uint32_t get_rounded_objective_value() {
        const double value = get_objective_value();
        assert(value >= 0);
        return static_cast<uint32_t>(lround(value));
    }

    /// @brief returns the objective coefficient of column j
    inline double get_objective_coefficient(const std::size_t j) {
        assert(j < get_nof_cols());
        return lp.getLp().col_cost_[j];
    }

    /**
     * @brief returns the lower and upper bound offset for 3-cycle inequalities ijk
     * due to permanently fixed variables that are not incorporated in the lp
     */
    inline int get_3cycle_bound_offset(const std::size_t ij, const std::size_t jk, const std::size_t ik) {
        assert(ij < n_choose_2);
        assert(jk < n_choose_2);
        assert(ik < n_choose_2);
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
    inline std::pair<double, double> get_row_bounds(const std::size_t i) {
        assert(i < get_nof_rows());
        const HighsLp &lp_ = lp.getLp();
        const double lb = lp_.row_lower_[i];
        const double ub = lp_.row_upper_[i];
        return std::make_pair(lb, ub);
    }

    /// @brief returns the value of row i
    inline double get_row_value(const std::size_t i) { return lp.getSolution().row_value[i]; }

    //
    // information methods
    //

    /**
     * @brief returns has_row_lb(i) && x > 0.
     *
     * @param i row index
     * @return has_row_lb(i) && x > 0.
     */
    inline bool has_row_slack(const std::size_t i) {
        const double x = get_row_value(i);
        const auto [lb, ub] = get_row_bounds(i);
        return x > lb + PACE_CONST_FEASIBILITY_TOLERANCE && x < ub - PACE_CONST_FEASIBILITY_TOLERANCE;
    }

    /**
     * @brief checks if a column variable of the lp is integral
     *
     * @param j column index
     * @return true if integral (that is, it is in the params.tol_bnd open neighborhood of an
     * integer)
     * @return false otherwise
     */
    inline bool is_column_integral(const std::size_t j) {
        assert(j < get_nof_cols());
        const double x = lp.getSolution().col_value[j];
        constexpr double ub = 1. - PACE_CONST_INTEGER_TOLERANCE;
        if (x > PACE_CONST_INTEGER_TOLERANCE && x < ub) {
            return false;
        }
        return true;
    }

    /**
     * @brief returns a value based on the integrality of the current solution
     *
     * @return true if integral
     * @return false otherwise
     */
    bool is_integral() {
        for (std::size_t j = 0; j < get_nof_cols(); ++j) {
            if (!is_column_integral(j)) {
                return false;
            }
        }
        return true;
    }

    /// @brief returns if an optimal feasible solution has been found
    bool is_optimal() {
        const HighsModelStatus &model_status = lp.getModelStatus();
        return model_status == HighsModelStatus::kOptimal;
    }

    /// @brief returns value of x < -1e-7
    inline bool is_3cycle_lb_violated(const double &x) { return x < -PACE_CONST_FEASIBILITY_TOLERANCE; }

    /// @brief returns value of x > 1 + 1e-7
    inline bool is_3cycle_ub_violated(const double &x) {
        constexpr double ub = 1. + PACE_CONST_FEASIBILITY_TOLERANCE;
        return x > ub;
    }

    //
    // 3-cycle methods
    //

    /**
     * @brief expects the violated 3-cycle ieqs in buckets.
     * adds the most violated <= nof_max_new_rows to the lp.
     *
     * @return std::size_t number of new rows
     */
    inline std::size_t add_3cycle_rows() {
        const std::size_t nof_new_rows = std::min(get_nof_bucket_entries(), nof_max_new_rows);
        if (nof_new_rows <= 0) return 0;

        lower_bounds.clear();
        upper_bounds.clear();
        starts.clear();
        indices.clear();
        values.clear();

        std::size_t i = 0;
        bool go_on = true;
        for (auto r_it = buckets.rbegin(); r_it != buckets.rend() && go_on; ++r_it) {
            for (const auto &[u, v, w] : *r_it) {
                const std::size_t uv = flat_index(n, n_choose_2, u, v);
                const std::size_t vw = flat_index(n, n_choose_2, v, w);
                const std::size_t uw = flat_index(n, n_choose_2, u, w);

                const double bound_offset = get_3cycle_bound_offset(uv, vw, uw);
                lower_bounds.emplace_back(0. + bound_offset);
                upper_bounds.emplace_back(1. + bound_offset);
                starts.emplace_back(indices.size());
                if (magic[uv] >= 0) {
                    indices.emplace_back(magic[uv]);
                    values.emplace_back(1.);
                }
                if (magic[uw] >= 0) {
                    indices.emplace_back(magic[uw]);
                    values.emplace_back(-1.);
                }
                if (magic[vw] >= 0) {
                    indices.emplace_back(magic[vw]);
                    values.emplace_back(1.);
                }
                assert(magic[uv] >= 0 || magic[uw] >= 0 || magic[vw] >= 0);

                rows.emplace_back(u, v, w);

                ++i;
                if (i >= nof_max_new_rows) {
                    go_on = false;
                    break;
                }
            }
        }
        assert(i == nof_new_rows);

        info.nof_added_rows = nof_new_rows;
        lp.addRows(nof_new_rows, &lower_bounds[0], &upper_bounds[0],  //
                   indices.size(), &starts[0], &indices[0], &values[0]);
        return nof_new_rows;
    }

    /**
     * @brief checks if the 3-cycle inequality for ijk/ikj is violated
     *
     * @return true if so
     * @return false otherwise
     */
    inline bool check_3cycle(const T &u, const T &v, const T &w) {
        const std::size_t uv = flat_index(n, n_choose_2, u, v);
        const std::size_t vw = flat_index(n, n_choose_2, v, w);
        const std::size_t uw = flat_index(n, n_choose_2, u, w);
        assert(uv < uw);
        assert(uw < vw);

        const double x = get_3cycle_value(uv, vw, uw);
        constexpr double interval_width = 1. + 2e-7;
        if (is_3cycle_lb_violated(x)) {
            const double x_normalized = -x / interval_width;
            get_bucket(x_normalized).emplace_back(u, v, w);
            return true;
        }
        if (is_3cycle_ub_violated(x)) {
            constexpr double ub = 1. + PACE_CONST_FEASIBILITY_TOLERANCE;
            const double x_normalized = (x - ub) / interval_width;
            get_bucket(x_normalized).emplace_back(u, v, w);
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
        clear_buckets();
        bool stop = false;
        // yes, I use the dark forces here, because it speeds things up and is imo cleaner
        T u = u_old, v = v_old, w = w_old;
        goto check_3cycles_in_for;
        for (; u < n - 2; ++u) {
            for (v = u + 1; v < n - 1; ++v) {
                for (w = v + 1; w < n; ++w) {
                check_3cycles_in_for:  // to pick off where we left
                    assert(u < v);
                    assert(v < w);
                    check_3cycle(u, v, w);
                    stop = is_last_bucket_full();
                    if (stop) goto check_3cycles_after_for;
                }
            }
        }
        ++info.nof_iterations_3cycles;
        if (nof_max_new_rows < 4 * PACE_CONST_NOF_CYCLE_CONSTRAINTS) {
            nof_max_new_rows *= 2;
        }
        for (u = 0; u < n - 2; ++u) {
            bool u_eq_old = u == u_old;
            for (v = u + 1; v < n - 1; ++v) {
                bool v_eq_old = v == v_old;
                for (w = v + 1; w < n; ++w) {
                    if (u_eq_old && v_eq_old && w == w_old) goto check_3cycles_after_for;
                    assert(u < v);
                    assert(v < w);
                    check_3cycle(u, v, w);
                    stop = is_last_bucket_full();
                    if (stop) goto check_3cycles_after_for;
                }
            }
        }

    check_3cycles_after_for:
        info.nof_bucket_entries = get_nof_bucket_entries();
        u_old = u, v_old = v, w_old = w;
        return add_3cycle_rows() > 0;
    }

    /**
     *  @brief returns x_uv + x_vw - x_uw of the lp
     */
    inline double get_3cycle_value(const std::size_t &uv, const std::size_t &vw, const std::size_t &uw) {
        const double x_uv = get_variable_value(uv);
        const double x_vw = get_variable_value(vw);
        const double x_uw = get_variable_value(uw);
        return x_uv + x_vw - x_uw;
    }

    //
    // bucket methods
    //

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
    inline std::vector<bucket_entry> &get_bucket(const double val) {
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
     * @return true if last bucket has >= nof_max_new_rows elements
     * @return false otherwise
     */
    inline bool is_last_bucket_full() { return (*buckets.rbegin()).size() >= nof_max_new_rows; }
};

};  // namespace pace

#endif
