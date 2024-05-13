#ifndef PACE_EXACT_HIGHS_BASE_HPP
#define PACE_EXACT_HIGHS_BASE_HPP

#include "Highs.h"
#include "model/instance.hpp"

namespace pace {

class highs_base : public instance_view {
   protected:
    Highs solver;
    HighsStatus status{HighsStatus::kOk};

    /*
     * for row adding
     */
    std::vector<double> lower_bounds;
    std::vector<double> upper_bounds;
    std::vector<HighsInt> starts;
    std::vector<HighsInt> indices;
    std::vector<double> values;

    highs_base(instance &instance_) : instance_view(instance_) {
        status = solver.setOptionValue("presolve", "off");
        assert(status == HighsStatus::kOk);
        status = solver.setOptionValue("solver", "simplex");
        assert(status == HighsStatus::kOk);
        status = solver.setOptionValue("parallel", "off");
        assert(status == HighsStatus::kOk);
        status = solver.setOptionValue("log_to_console", false);
        assert(status == HighsStatus::kOk);

        PACE_DEBUG_PRINTF("start add_columns\n");
        add_columns();
        PACE_DEBUG_PRINTF("end   add_columns\n");
    }

    /*
     * column methods
     */

    /**
     * @brief adds all variables of the ilp formulation of one-sided crossing minimization
     * and sets their objective coefficients from the crossing numbers
     */
    inline void add_columns() {
        crossing_number_t obj_offset = objective_offset();

        // add variables which are not yet settled
        const crossing_matrix &cr{cr_matrix()};
        for (const auto &[u, v] : unsettled_pairs()) {
            const HighsInt l = solver.getNumCol();
            assert(l == magic()[flat_index(n_free, n_free_2, u, v)]);
            add_column();

            const crossing_number_t &c_uv = cr(u, v);
            const crossing_number_t &c_vu = cr(v, u);
            obj_offset += c_vu;
            if (c_uv != c_vu) {
                change_column_cost(l, static_cast<double>(c_uv) - static_cast<double>(c_vu));
            }
        }

        // set constant term (shift/offset) in the objective function
        solver.changeObjectiveOffset(obj_offset);
        assert(get_n_cols() <= n_free_2);
    }

    /**
     * @brief adds a column with bounds 0. <= (*) <= 1.
     */
    inline void add_column() {
        status = solver.addVar(0., 1.);
        assert(status == HighsStatus::kOk);
    }

    /**
     * @brief changes cost of column j
     */
    inline void change_column_cost(const std::size_t j, const double cost) {
        assert(j < get_n_cols());
        status = solver.changeColCost(j, cost);
        assert(status == HighsStatus::kOk);
    }

    /*
     * row methods
     */

    /**
     * @brief prepares and adds entries to auxiliary vectors for adding 3-cycle ieq uvw/uwv to lp
     */
    inline void add_3cycle_row_to_aux_vectors(const vertex_t &u, const vertex_t &v, const vertex_t &w) {
        // at least two must be in the lp as variables since we compute transitive hull in oracle
        assert(get_n_vars_in_lp(u, v, w) >= 2);
        assert(u < v);
        assert(v < w);

        const auto [uv, vw, uw] = flat_indices(n_free, n_free_2, u, v, w);
        const double bound_offset = get_3cycle_bound_offset(uv, vw, uw);
        const std::vector<magic_t> &m = magic();

        lower_bounds.emplace_back(0. + bound_offset);
        upper_bounds.emplace_back(1. + bound_offset);
        starts.emplace_back(indices.size());
        if (m[uv] >= 0) {
            indices.emplace_back(m[uv]);
            values.emplace_back(1.);
        }
        if (m[uw] >= 0) {
            indices.emplace_back(m[uw]);
            values.emplace_back(-1.);
        }
        if (m[vw] >= 0) {
            indices.emplace_back(m[vw]);
            values.emplace_back(1.);
        }
    }

    /*
     * getter
     */

   public:
    void copy_column_values(std::vector<double> &col_value) const { col_value = solver.getSolution().col_value; }

    const std::vector<double> &get_column_values() { return solver.getSolution().col_value; }

    /**
     * @brief returns value of column j
     */
    double get_column_value(const std::size_t j) {
        assert(j < get_n_cols());
        return solver.getSolution().col_value[j];
    }

    /**
     * @brief returns value of variable x_uv
     */
    double get_variable_value(const std::size_t uv) {
        assert(uv < n_free_2);
        const std::vector<magic_t> &m{magic()};
        if (m[uv] < 0) {
            return static_cast<double>(~m[uv]);
        } else {
            return solver.getSolution().col_value[m[uv]];
        }
    }

    /**
     * @brief returns the number of rows
     */
    std::size_t get_n_cols() const { return solver.getNumCol(); }

    /**
     * @brief returns the number of rows
     */
    std::size_t get_n_rows() const { return solver.getNumRow(); }

    /**
     * @brief returns the objective value of the lp
     */
    double get_objective_value() const { return solver.getInfo().objective_function_value; }

    /**
     * @brief returns the objective value of the lp rounded to the next integer
     */
    crossing_number_t get_rounded_objective_value() const {
        const double value = get_objective_value();
        assert(value >= 0);
        return static_cast<crossing_number_t>(lround(value));
    }

    /**
     * @brief returns the objective coefficient (cost) of column j
     */
    inline double get_objective_coefficient(const std::size_t j) {
        assert(j < get_n_cols());
        return solver.getLp().col_cost_[j];
    }

    /**
     * @brief returns the lower and upper bound offset for 3-cycle inequalities ijk
     * due to permanently fixed variables that are not incorporated in the lp
     */
    inline int  //
    get_3cycle_bound_offset(const std::size_t uv, const std::size_t vw, const std::size_t uw) {
        assert(uv < n_free_2);
        assert(vw < n_free_2);
        assert(uw < n_free_2);
        double offset = 0.;
        const std::vector<magic_t> &m{magic()};
        if (m[uv] < 0) offset -= ~m[uv];
        if (m[uw] < 0) offset += ~m[uw];
        if (m[vw] < 0) offset -= ~m[vw];
        return offset;
    }

    /**
     * @brief returns the lower and upper bound offset for 3-cycle inequality with index i
     * due to permanently fixed variables that are not incorporated in the lp
     */
    inline std::pair<double, double> get_row_bounds(const std::size_t i) {
        assert(i < get_n_rows());
        const HighsLp &lp_ = solver.getLp();
        const double lb = lp_.row_lower_[i];
        const double ub = lp_.row_upper_[i];
        return std::make_pair(lb, ub);
    }

    /**
     * @brief returns the value of row i
     */
    inline double get_row_value(const std::size_t i) {
        assert(i < get_n_rows());
        return solver.getSolution().row_value[i];
    }

    /**
     * @brief returns number of variables x_uv, x_vw, x_uw are in the lp
     */
    inline int8_t get_n_vars_in_lp(const vertex_t &u, const vertex_t &v, const vertex_t &w) {
        const auto [uv, vw, uw] = flat_indices(n_free, n_free_2, u, v, w);
        uint8_t n_vars_in_lp{0};
        const std::vector<magic_t> m{magic()};
        n_vars_in_lp += m[uv] >= 0;
        n_vars_in_lp += m[uw] >= 0;
        n_vars_in_lp += m[vw] >= 0;
        return n_vars_in_lp;
    }

    /*
     * boolean
     */

    /**
     * @brief returns true iff an optimal feasible solution has been found
     */
    bool is_optimal() {
        const HighsModelStatus &model_status = solver.getModelStatus();
        return model_status == HighsModelStatus::kOptimal;
    }

    /**
     * @brief returns true iff solution is feasible
     */
    bool is_feasible() {
        const HighsModelStatus &model_status = solver.getModelStatus();
        return model_status != HighsModelStatus::kInfeasible;
    }

    /**
     * @brief returns true iff solutions is infeasible
     */
    bool is_infeasible() {
        const HighsModelStatus &model_status = solver.getModelStatus();
        return model_status == HighsModelStatus::kInfeasible;
    }
};

};  // namespace pace

#endif
