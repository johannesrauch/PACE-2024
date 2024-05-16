#include "exact/highs_base.hpp"

namespace pace {

//
// ctor
//

highs_base::highs_base(instance &instance_) : instance_view(instance_) {
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

//
// column methods
//

void highs_base::add_columns() {
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
            change_column_cost(
                l, static_cast<double>(c_uv) - static_cast<double>(c_vu));
        }
    }

    // set constant term (shift/offset) in the objective function
    solver.changeObjectiveOffset(obj_offset);
    assert(get_n_cols() <= n_free_2);
}

//
// row methods
//

void highs_base::add_3cycle_row_to_aux_vectors(  //
    const vertex_t &u, const vertex_t &v, const vertex_t &w) {
    // at least two must be in the lp as variables since we compute transitive
    // hull in oracle
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

};  // namespace pace
