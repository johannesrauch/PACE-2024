#ifndef PACE_EXACT_LP_WRAPPER_HPP
#define PACE_EXACT_LP_WRAPPER_HPP

#include <unordered_map>
#include <unordered_set>

#include "exact/highs_base.hpp"
#include "exact/info_structs.hpp"
#include "utils/index_utils.hpp"
#include "utils/randomness_utils.hpp"

namespace pace {

struct highs_lp_params {
    std::size_t max_new_rows{
        256};  ///< maximum number of new rows per check_3cycles call
    const uint16_t max_initial_rows{16384};  ///< maximum number of initial rows
    const uint8_t max_delete_rows_3cycle_iters{
        64};  ///< maximum number of 3-cycle iterations with row deletion
    const uint16_t max_initial_solve_iters{200};
    const uint8_t max_initial_solve_3cycle_iters{
        2};  ///< maximum number of times max_new_rows is doubled
    const int32_t max_initial_solve_simplex_iters{30000};

    const uint8_t n_buckets{8};

    const double tol_feasibility{1e-7};
    const double tol_integer{1e-6};
};

/**
 * @brief wrapper class for highs solver; creates and manages ilp relaxation of
 * the instance
 */
class highs_lp : public highs_base {
    highs_lp_info info;
    highs_lp_params params;

    /**
     * @brief number of rows before cut() added new rows
     */
    std::size_t n_old_rows{0};

    /**
     * @brief stores indices of the to deleted rows in
     * delete_positive_slack_rows()
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
    std::vector<std::vector<triple>> buckets{params.n_buckets};

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
            return hash(std::get<0>(t)) ^ hash(std::get<1>(t)) ^
                   hash(std::get<2>(t));
        }
    };

    std::unordered_map<triple, row_info, triple_hash> rows_info;

    std::list<triple> rows;

   public:
    highs_lp(instance &instance_, highs_lp_params params = highs_lp_params());

    highs_lp(const highs_lp &rhs) = delete;
    highs_lp(highs_lp &&rhs) = delete;
    highs_lp &operator=(const highs_lp &rhs) = delete;
    highs_lp &operator=(highs_lp &&rhs) = delete;

    //
    // lp solving methods
    //

    /**
     * @brief search for violated constraints and add these to the lp
     *
     * @return true if violated constraints found
     * @return false otherwise
     */
    bool cut();

    inline void set_simplex_iteration_limit(const int32_t max_simplex_iter) {
        status =
            solver.setOptionValue("simplex_iteration_limit", max_simplex_iter);
        assert(status == HighsStatus::kOk);
    }

    /**
     * @brief solves the current version of the lp relaxation
     * with an upper limit of the number of simplex iterations
     *
     * @param max_simplex_iter max simplex iterations
     */
    void run(
        const int32_t max_simplex_iter = std::numeric_limits<int32_t>::max(),
        const bool bookkeeping = true);

    void initial_partial_solve();

    void resolve();

    //
    // row modification methods
    //

    /**
     * @brief delete rows with positive slack from the lp
     */
    void delete_positive_slack_rows();

    //
    // column modification methods
    //

    /**
     * @brief changes upper and lower bound of column j
     */
    inline void change_column_bounds(const std::size_t j, const double lb,
                                     const double ub) {
        assert(j < get_n_cols());
        status = solver.changeColBounds(j, lb, ub);
        assert(status == HighsStatus::kOk);
    }

    /**
     * @brief fixes column j to fix_to
     */
    inline void fix_column(const std::size_t j, const double fix_to) {
        assert(0. <= fix_to);
        assert(fix_to <= 1.);
        PACE_DEBUG_PRINTF("\tfixed column %5d to %1.0f\n", j, fix_to);
        change_column_bounds(j, fix_to, fix_to);
    }

    /**
     * @brief fixes columns based on the following:
     * if lower_bound = sum min(c_ij, c_ji), where c_ij and c_ji are crossing
     * numbers, and coeff = c_ij - c_ji is the objective coefficient of variable
     * x_ij, and abs(c_ij - c_ji) >= upper_bound - lower_bound, then we are able
     * to fix x_ij to 0 if coeff > 0 and to 1 otherwise
     */
    void fix_columns();

    /**
     * @brief resets bounds of column j to 0 <= . <= 1
     */
    inline void unfix_column(const std::size_t j) {
        PACE_DEBUG_PRINTF("\tunfixed column %5d\n", j);
        change_column_bounds(j, 0., 1.);
    }

    //
    // getter
    //

    const highs_lp_info &get_info() { return info; }

    //
    // information methods
    //

    /**
     * @brief returns true iff row i is in the open interval (0; 1)
     */
    inline bool has_row_slack(const std::size_t i) {
        const double x = get_row_value(i);
        const auto [lb, ub] = get_row_bounds(i);
        return x > lb + params.tol_feasibility &&
               x < ub - params.tol_feasibility;
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

    inline bool is_variable_integral(const std::size_t uv) {
        assert(uv < n_free_2);
        const magic_t j = magic()[uv];
        if (j < 0) return true;
        return is_column_integral(j);
    }

    /**
     * @brief returns true iff current solution is integral
     */
    bool is_integral();

    /**
     * @brief returns x < -params.tol_feasibility
     */
    inline bool is_3cycle_lb_violated(const double &x) {
        return x < -params.tol_feasibility;
    }

    /**
     * @brief returns x > 1 + params.tol_feasibility
     */
    inline bool is_3cycle_ub_violated(const double &x) {
        return x > 1. + params.tol_feasibility;
    }

    /**
     * @brief returns true iff 3-cycle ieq for u < v < w is interesting,
     * that is, the number of crossings of a combination that violates a 3-cycle
     * ieq is less than the number of crossings of a combination (a permutation
     * of uvw) that does not
     */
    bool is_3cycle_interesting(const vertex_t &u, const vertex_t &v,
                               const vertex_t &w);

    //
    // 3-cycle methods
    //

    bool add_3cycle_row_to_internal_rows(  //
        const vertex_t &u, const vertex_t &v, const vertex_t &w);

    void add_initial_rows_from_ordering();

    /**
     * @brief adds at most params.max_initial_rows "interesting" rows to the lp
     */
    void add_initial_rows_from_interesting();

    /**
     * @brief expects the violated 3-cycle ieqs in buckets.
     * adds the most violated <= params.max_new_rows to the lp.
     *
     * @return std::size_t number of new rows
     */
    std::size_t add_3cycle_rows();

    /**
     * @brief checks if the 3-cycle inequalities of u < v < w is violated,
     * adds them to buckets in positive case
     *
     * @return true if so
     * @return false otherwise
     */
    bool check_3cycle(const vertex_t &u, const vertex_t &v, const vertex_t &w);

    /**
     * @brief checks if the 3-cycle inequalities are violated
     *
     * @return true if so
     * @return false otherwise
     */
    bool check_3cycles();
    bool check_3cycles_depr();

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
    inline double get_3cycle_value(const std::size_t &uv, const std::size_t &vw,
                                   const std::size_t &uw) {
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
        info.min_viol_score = std::numeric_limits<double>::max();
        info.max_viol_score = std::numeric_limits<double>::min();
    }

    /**
     * @brief given by how much the 3-cycle inequality is violated,
     * returns the corresponding bucket
     *
     * @param val in [0,1)
     * @return std::vector<triple>& the corresponding bucket
     */
    inline std::vector<triple> &get_bucket(double val) {
        assert(0 <= val);
        // assert(val < 1);
        val = std::min(val, 0.99);
        const std::size_t i = static_cast<std::size_t>(val * params.n_buckets);
        assert(i < params.n_buckets);
        return buckets[i];
    }

    inline double get_objective_value_contribution(const vertex_t &u,
                                                   const vertex_t &v) {
        const std::size_t uv = flat_index(n_free, n_free_2, u, v);
        const crossing_matrix &cr = cr_matrix();
        const std::vector<magic_t> &m = magic();
        return (cr(u, v) - cr(v, u)) * get_variable_value(uv);
    }

    inline std::vector<triple> &get_bucket_from_score(const vertex_t &u,
                                                      const vertex_t &v,
                                                      const vertex_t &w) {
        const double score = get_objective_value_contribution(u, v) -  //
                             get_objective_value_contribution(u, w) +  //
                             get_objective_value_contribution(v, w);
        info.min_viol_score = std::min(score, info.min_viol_score);
        info.max_viol_score = std::max(score, info.max_viol_score);
        const double interval = info.max_viol_score - info.min_viol_score;
        double val = 0.99;
        if (interval >= 0.1) val -= (score - info.min_viol_score) / interval;
        val = std::max(0., val);
        return get_bucket(val);
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
     * @return true if last bucket has >= params.max_new_rows elements
     * @return false otherwise
     */
    inline bool is_last_bucket_full() {
        return (*buckets.rbegin()).size() >= params.max_new_rows;
    }

    inline bool is_n_bucket_entries_large() {
        return get_n_bucket_entries() >= 10 * params.max_new_rows;
    }

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
    void update_3cycle_iteration_info();

    void update_simplex_info();

    inline void reset_row_info() {
        info.n_bucket_entries = 0;
        info.n_added_rows = 0;
        info.n_deleted_rows = 0;
        info.n_delete_rows_spared = 0;
        info.tried_deleting_rows = false;
    }
};

};  // namespace pace

#endif
