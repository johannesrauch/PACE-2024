#ifndef PACE_EXACT_LP_WRAPPER_HPP
#define PACE_EXACT_LP_WRAPPER_HPP

#include <unordered_map>

#include "exact/highs_base.hpp"
#include "exact/info_structs.hpp"
#include "utils/index_utils.hpp"
#include "utils/randomness_utils.hpp"

#ifndef PACE_CONST_N_MAX_NEW_ROWS
#define PACE_CONST_N_MAX_NEW_ROWS 256u
#endif

#ifndef PACE_CONST_N_MAX_INIT_ROWS
#define PACE_CONST_N_MAX_INIT_ROWS 8192
#endif

#ifndef PACE_CONST_N_BUCKETS
#define PACE_CONST_N_BUCKETS 8u
#endif

#ifndef PACE_CONST_TOL_FEASIBILITY
#define PACE_CONST_TOL_FEASIBILITY 1e-7
#endif

#ifndef PACE_CONST_TOL_INTEGER
#define PACE_CONST_TOL_INTEGER 1e-6
#endif

namespace pace {

struct highs_wrapper_params {
    std::size_t limit_new_rows{PACE_CONST_N_MAX_NEW_ROWS};
    const uint8_t limit_new_rows_double{2};
    const std::size_t limit_initial_rows{PACE_CONST_N_MAX_INIT_ROWS};
    const uint8_t limit_delete_slack_row{16};
    const uint8_t limit_delete_rows{64};

    const uint8_t n_buckets{PACE_CONST_N_BUCKETS};

    const double tol_feasibility{PACE_CONST_TOL_FEASIBILITY};
    const double tol_integer{PACE_CONST_TOL_INTEGER};

    double p_delete_slack_row{0.};  // should be < 0.001 if positive
};

/**
 * @brief wrapper class for highs solver; creates and manages ilp relaxation of the instance
 */
class highs_lp : public highs_base {
    highs_wrapper_info info;
    highs_wrapper_params params;

    /**
     * @brief number of rows before cut() added new rows
     */
    std::size_t n_old_rows{0};

    /**
     * @brief the first n_fix_rows are not subject to deletion
     */
    std::size_t n_fix_rows{0};

    /**
     * @brief stores indices of the to deleted rows in delete_positive_slack_rows()
     */
    std::vector<HighsInt> rows_to_delete;

    /**
     * @brief where we last stopped checking 3-cycle ieqs
     */
    vertex_t u_old{0}, v_old{1}, w_old{2};

    //
    // bucket attributes and methods
    //

    using triple = std::tuple<vertex_t, vertex_t, vertex_t>;

    /**
     * @brief buckets for bucket sorting violated 3-cycle inequalities
     */
    std::vector<std::vector<triple>> buckets{PACE_CONST_N_BUCKETS};

    //
    // for row management, todo: unify bucket_entry and triple
    //

    struct row_info {
        uint16_t n_deleted{0};
        uint16_t n_spared{0};
    };

    struct triple_hash {
        std::size_t operator()(const triple &t) const {
            const std::hash<vertex_t> hash;
            return hash(std::get<0>(t)) ^ hash(std::get<1>(t)) ^ hash(std::get<2>(t));
        }
    };

    std::unordered_map<triple, row_info, triple_hash> rows_info;

    std::list<triple> rows;

   public:
    highs_lp(instance &instance_, highs_wrapper_params params = highs_wrapper_params())
        : highs_base(instance_), info{u_old, v_old, w_old}, params(params) {
        rows_to_delete.reserve(params.limit_new_rows);
        lower_bounds.reserve(params.limit_new_rows);
        upper_bounds.reserve(params.limit_new_rows);
        starts.reserve(params.limit_new_rows);
        indices.reserve(3 * params.limit_new_rows);
        values.reserve(3 * params.limit_new_rows);

        info.n_cols = get_n_cols();
        // PACE_DEBUG_PRINTF("start add_initial_rows\n");
        // add_initial_rows();
        // PACE_DEBUG_PRINTF("end   add_initial_rows\n");
    }

    highs_lp(const highs_lp &rhs) = delete;
    highs_lp(highs_lp &&rhs) = delete;
    highs_lp &operator=(const highs_lp &rhs) = delete;
    highs_lp &operator=(highs_lp &&rhs) = delete;

    ~highs_lp() {}

    //
    // lp solving methods
    //

    /**
     * @brief search for violated constraints and add these to the lp
     *
     * @return true if violated constraints found
     * @return false otherwise
     */
    bool cut() {
        n_old_rows = get_n_rows();
        bool success = check_3cycles();
        return success;
    }

    inline void set_simplex_iteration_limit(const int32_t limit_simplex_it) {
        status = solver.setOptionValue("simplex_iteration_limit", limit_simplex_it);
        assert(status == HighsStatus::kOk);
    }

    /**
     * @brief solves the current version of the lp relaxation
     */
    void run() {
        solver.run();
        update_simplex_info();
    }

    /**
     * @brief for reliability branching; no bookkeeping
     *
     * @param limit_simplex_it max simplex iterations
     */
    void run(const int32_t limit_simplex_it) {
        set_simplex_iteration_limit(limit_simplex_it);
        solver.run();
        set_simplex_iteration_limit(std::numeric_limits<int32_t>::max());
    }

    //
    // row modification methods
    //

    /**
     * @brief delete rows with positive slack from the lp
     */
    void delete_positive_slack_rows() {
        reset_delete_count_info();
        if (info.n_iterations_3cycles > params.limit_delete_rows) return;

        // gather rows to delete
        rows_to_delete.clear();
        auto it_rows = rows.begin();
        for (std::size_t i = n_fix_rows; i < n_old_rows; ++i) {
            assert(it_rows != rows.end());
            auto it = rows_info.find(*it_rows);
            assert(it != rows_info.end());

            const bool has_slack = has_row_slack(i);
            const bool spare = it->second.n_deleted > it->second.n_spared;
            const bool remove = has_slack && !spare;
            const bool bad_luck =
                it->second.n_deleted < params.limit_delete_slack_row && coinflip(params.p_delete_slack_row);

            if (remove || bad_luck) {
                rows_to_delete.emplace_back(i);
                ++it->second.n_deleted;
                it->second.n_spared = 0;
                it_rows = rows.erase(it_rows);
            } else {
                ++it_rows;
            }

            if (has_slack && spare) {
                ++it->second.n_spared;
                ++info.n_delete_rows_spared;
            }
            info.n_deleted_rows_slack += remove;
            info.n_deleted_rows_bad_luck += !remove && bad_luck;
        }

        info.n_deleted_rows = rows_to_delete.size();
        if (info.n_deleted_rows > 0) {
            solver.deleteRows(info.n_deleted_rows, &rows_to_delete[0]);
        }
    }

    //
    // column modification methods
    //

    /**
     * @brief changes upper and lower bound of column j
     */
    inline void change_column_bounds(const std::size_t j, const double lb, const double ub) {
        assert(j < get_n_cols());
        status = solver.changeColBounds(j, lb, ub);
        assert(status == HighsStatus::kOk);
    }

    /**
     * @brief fixes column j to fix_to
     */
    void fix_column(const std::size_t j, const double fix_to) {
        assert(0. <= fix_to);
        assert(fix_to <= 1.);
        PACE_DEBUG_PRINTF("\tfixed column %5d to %1.0f\n", j, fix_to);
        change_column_bounds(j, fix_to, fix_to);
    }

    /**
     * @brief fixes columns based on the following:
     * if lower_bound = sum min(c_ij, c_ji), where c_ij and c_ji are crossing numbers,
     * and coeff = c_ij - c_ji is the objective coefficient of variable x_ij,
     * and abs(c_ij - c_ji) >= upper_bound - lower_bound,
     * then we are able to fix x_ij to 0 if coeff > 0 and to 1 otherwise
     */
    void fix_columns() {
        for (std::size_t j = 0; j < get_n_cols(); ++j) {
            const crossing_number_t diff = upper_bound - lower_bound();
            const double coeff = get_objective_coefficient(j);

            if (std::abs(coeff) >= static_cast<double>(diff) && coeff != 0.) {
                fix_column(j, coeff > 0 ? 0. : 1.);
                PACE_DEBUG_PRINTF("(optimality condition)\n");
            }
        }
    }

    /**
     * @brief resets bounds of column j to 0 <= . <= 1
     */
    void unfix_column(const std::size_t j) {
        PACE_DEBUG_PRINTF("\tunfixed column %5d\n", j);
        change_column_bounds(j, 0., 1.);
    }

    //
    // getter
    //

    const highs_wrapper_info &get_info() { return info; }

    //
    // information methods
    //

    /**
     * @brief returns true iff row i is in the open interval (0; 1)
     */
    inline bool has_row_slack(const std::size_t i) {
        const double x = get_row_value(i);
        const auto [lb, ub] = get_row_bounds(i);
        return x > lb + params.tol_feasibility && x < ub - params.tol_feasibility;
    }

    /**
     * @brief returns true iff a column j of the lp is integral
     */
    inline bool is_column_integral(const std::size_t j) {
        const double x = get_column_value(j);
        if (x > params.tol_integer && x < 1. - params.tol_integer) {
            return false;
        }
        return true;
    }

    /**
     * @brief returns true iff current solution is integral
     */
    bool is_integral() {
        for (std::size_t j = 0; j < get_n_cols(); ++j) {
            if (!is_column_integral(j)) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief returns x < -params.tol_feasibility
     */
    inline bool is_3cycle_lb_violated(const double &x) { return x < -params.tol_feasibility; }

    /**
     * @brief returns x > 1 + params.tol_feasibility
     */
    inline bool is_3cycle_ub_violated(const double &x) { return x > 1. + params.tol_feasibility; }

    /**
     * @brief returns true iff 3-cycle ieq for u < v < w is interesting,
     * that is, the number of crossings of a combination that violates a 3-cycle ieq
     * is less than the number of crossings of a combination (a permutation of uvw) that does not
     */
    inline bool is_3cycle_interesting(const vertex_t &u, const vertex_t &v, const vertex_t &w) {
        assert(u < v);
        assert(v < w);

        const crossing_matrix &cr = cr_matrix();

        const crossing_number_t c_uvw = cr(u, v) + cr(u, w) + cr(v, w);
        const crossing_number_t c_uwv = cr(u, v) + cr(u, w) + cr(w, v);
        const crossing_number_t c_vuw = cr(v, u) + cr(u, w) + cr(v, w);
        const crossing_number_t c_vwu = cr(v, u) + cr(w, u) + cr(v, w);
        const crossing_number_t c_wuv = cr(u, v) + cr(w, u) + cr(w, v);
        const crossing_number_t c_wvu = cr(v, u) + cr(w, u) + cr(w, v);
        const crossing_number_t c_min_allowed = std::min({c_uvw, c_uwv, c_vuw, c_vwu, c_wuv, c_wvu});

        const crossing_number_t c_uv_vw_wu = cr(u, v) + cr(v, w) + cr(w, u);
        const crossing_number_t c_vu_wv_uw = cr(v, u) + cr(w, v) + cr(u, w);
        const crossing_number_t c_min_forbidden = std::min(c_uv_vw_wu, c_vu_wv_uw);

        return c_min_forbidden < c_min_allowed;
    }

    //
    // 3-cycle methods
    //

    /**
     * @brief adds at most params.limit_initial_rows "interesting" rows to the lp
     */
    inline void add_initial_rows() {
        clear_aux_vectors();
        std::vector<triple> candidates;
        for (vertex_t u = 0u; u + 2u < n_free; ++u) {
            for (vertex_t v = u + 1u; v + 1u < n_free; ++v) {
                for (vertex_t w = v + 1u; w < n_free; ++w) {
                    assert(u < v);
                    assert(v < w);
                    if (get_n_vars_in_lp(u, v, w) >= 2 && is_3cycle_interesting(u, v, w)) {
                        candidates.emplace_back(u, v, w);
                    }
                }
            }
        }

        const double p = static_cast<double>(params.limit_initial_rows) / candidates.size();
        for (const auto &[u, v, w] : candidates) {
            if (coinflip(p)) {
                add_3cycle_row_to_aux_vectors(u, v, w);
                add_3cycle_row_to_internal_rows(u, v, w);
            }
        }

        solver.addRows(lower_bounds.size(), &lower_bounds[0], &upper_bounds[0],  //
                       indices.size(), &starts[0], &indices[0], &values[0]);

        info.n_init_rows_candidates = candidates.size();
        info.n_rows = get_n_rows();
    }

    inline void add_3cycle_row_to_internal_rows(const vertex_t &u, const vertex_t &v, const vertex_t &w) {
        rows.emplace_back(u, v, w);
        const triple uvw(u, v, w);
        auto it = rows_info.find(uvw);
        if (it == rows_info.end()) {
            rows_info.emplace(uvw, row_info());
        }
    }

    /**
     * @brief expects the violated 3-cycle ieqs in buckets.
     * adds the most violated <= params.limit_new_rows to the lp.
     *
     * @return std::size_t number of new rows
     */
    inline std::size_t add_3cycle_rows() {
        info.n_added_rows = std::min(get_n_bucket_entries(), params.limit_new_rows);
        if (info.n_added_rows <= 0) return 0;

        clear_aux_vectors();

        std::size_t i = 0;
        bool go_on = true;
        for (auto r_it = buckets.rbegin(); r_it != buckets.rend() && go_on; ++r_it) {
            for (const auto &[u, v, w] : *r_it) {
                add_3cycle_row_to_aux_vectors(u, v, w);
                add_3cycle_row_to_internal_rows(u, v, w);
                ++i;
                if (i >= params.limit_new_rows) {
                    go_on = false;
                    break;
                }
            }
        }
        assert(i == info.n_added_rows);

        solver.addRows(info.n_added_rows, &lower_bounds[0], &upper_bounds[0],  //
                       indices.size(), &starts[0], &indices[0], &values[0]);
        return info.n_added_rows;
    }

    /**
     * @brief checks if the 3-cycle inequalities of u < v < w is violated,
     * adds them to buckets in positive case
     *
     * @return true if so
     * @return false otherwise
     */
    inline bool check_3cycle(const vertex_t &u, const vertex_t &v, const vertex_t &w) {
        const auto [uv, vw, uw] = flat_indices(n_free, n_free_2, u, v, w);
        const double x = get_3cycle_value(uv, vw, uw);
        const double interval_width = 1. + 2. * params.tol_feasibility;
        if (is_3cycle_lb_violated(x)) {
            const double x_normalized = -x / interval_width;
            get_bucket(x_normalized).emplace_back(u, v, w);
            return true;
        }
        if (is_3cycle_ub_violated(x)) {
            const double ub = 1. + params.tol_feasibility;
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
        vertex_t u = u_old, v = v_old, w = w_old;
        goto check_3cycles_in_for;
        for (; u + 2u < n_free; ++u) {
            for (v = u + 1; v + 1u < n_free; ++v) {
                for (w = v + 1; w < n_free; ++w) {
                check_3cycles_in_for:  // to pick off where we left
                    assert(u < v);
                    assert(v < w);
                    check_3cycle(u, v, w);
                    stop = is_last_bucket_full();
                    if (stop) goto check_3cycles_after_for;
                }
            }
        }
        update_3cycle_iteration_info();
        for (u = 0; u + 2u < n_free; ++u) {
            bool u_eq_old = u == u_old;
            for (v = u + 1; v + 1u < n_free; ++v) {
                bool v_eq_old = v == v_old;
                for (w = v + 1; w < n_free; ++w) {
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
        info.n_bucket_entries = get_n_bucket_entries();
        u_old = u, v_old = v, w_old = w;
        return add_3cycle_rows() > 0;
    }

    /**
     * @brief clears auxiliary vectors for adding rows
     */
    inline void clear_aux_vectors() {
        lower_bounds.clear();
        upper_bounds.clear();
        starts.clear();
        indices.clear();
        values.clear();
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

    /**
     * @brief clears contents of each bucket in buckets
     */
    inline void clear_buckets() {
        for (std::vector<triple> &bucket : buckets) {
            bucket.clear();
        }
    }

    /**
     * @brief given by how much the 3-cycle inequality is violated,
     * returns the corresponding bucket
     *
     * @param val in (0,1]
     * @return std::vector<triple>& the corresponding bucket
     */
    inline std::vector<triple> &get_bucket(const double val) {
        assert(0 <= val);
        assert(val < 1);
        const std::size_t i = static_cast<std::size_t>(val * params.n_buckets);
        assert(i < params.n_buckets);
        return buckets[i];
    }

    /// @brief returns sum of the number of elements in each bucket of `buckets`
    inline std::size_t get_n_bucket_entries() {
        std::size_t n = 0;
        for (const auto &bucket : buckets) {
            n += bucket.size();
        }
        return n;
    }

    /**
     * @brief returns if the last bucket is full
     *
     * @return true if last bucket has >= params.limit_new_rows elements
     * @return false otherwise
     */
    inline bool is_last_bucket_full() { return (*buckets.rbegin()).size() >= params.limit_new_rows; }

    //
    // hot start methods
    //

    HighsInt freeze_basis() {
        HighsInt frozen_basis_id;
        status = solver.freezeBasis(frozen_basis_id);
        assert(status == HighsStatus::kOk);
        return frozen_basis_id;
    }

    void unfreeze_basis(const HighsInt frozen_basis_id) {
        status = solver.unfreezeBasis(frozen_basis_id);
        assert(status == HighsStatus::kOk);
    }

    //
    // bookkeeping methods
    //
   private:
    inline void reset_delete_count_info() {
        info.n_deleted_rows = 0;
        info.n_delete_rows_spared = 0;
        info.n_deleted_rows_bad_luck = 0;
        info.n_deleted_rows_slack = 0;
    }

    inline void update_3cycle_iteration_info() {
        if (info.n_iterations_3cycles < params.limit_new_rows_double) {
            params.limit_new_rows *= 2;
        }
        params.p_delete_slack_row /= 2;
        ++info.n_iterations_3cycles;
    }

    inline void update_simplex_info() {
        info.n_rows = get_n_rows();
        info.n_iterations_simplex = solver.getSimplexIterationCount();
        info.n_iterations_simplex_avg =
            (info.n_iterations_simplex_avg * info.n_runs + info.n_iterations_simplex) / (info.n_runs + 1);

        info.t_simplex = solver.getRunTime();

        info.objective_value = get_objective_value();
        ++info.n_runs;
    }
};

};  // namespace pace

#endif
