#ifndef PACE_EXACT_HIGHS_MIP_HPP
#define PACE_EXACT_HIGHS_MIP_HPP

#include "exact/highs_base.hpp"
#include "utils/crossings_utils.hpp"
#include "utils/restr_graph_utils.hpp"

namespace pace {

class highs_mip : public highs_base {
   public:
    highs_mip(instance &instance_) : highs_base(instance_) {
        set_columns_integral();
        add_all_3cycle_rows();
    }

    /**
     * @brief solves the oscm instance exactly with highs branch and cut mip solver
     *
     * @param ordering out parameter
     * @return crossing_number_t number of crossings
     */
    crossing_number_t operator()(std::vector<vertex_t> &ordering) {
        PACE_DEBUG_PRINTF("\n");
        ordering.resize(n_free);
        for (vertex_t i = 0; i < n_free; ++i) ordering[i] = i;

        crossing_number_t n_crossings;
        if (get_n_cols() > 0) {
            solver.run();
            assert(is_optimal());
            n_crossings = get_rounded_objective_value();
#ifndef NDEBUG
            const bool acyclic =
#endif
                build_restr_graph_ordering(  //
                    solver.getSolution().col_value, unsettled_pairs(), 1e-6, restriction_graph(), ordering);
            assert(acyclic);
        } else {
            topological_sort(restriction_graph(), ordering);
            n_crossings = number_of_crossings(graph, ordering);
        }
        update_ordering(ordering, n_crossings);

        return n_crossings;
    }

   private:
    /**
     * @brief adds all relevant 3-cycle ieq to lp
     */
    inline void add_all_3cycle_rows() {
        std::size_t n_free_3 = n_free_2 * (n_free / 3);  // approx n_free choose 3
        lower_bounds.reserve(n_free_3);
        upper_bounds.reserve(n_free_3);
        starts.reserve(n_free_3);
        indices.reserve(n_free_3);
        values.reserve(n_free_3);

        for (vertex_t u = 0; u + 2u < n_free; ++u) {
            for (vertex_t v = u + 1u; v + 1u < n_free; ++v) {
                for (vertex_t w = v + 1u; w < n_free; ++w) {
                    assert(u < v);
                    assert(v < w);
                    if (get_n_vars_in_lp(u, v, w) >= 2) {
                        add_3cycle_row_to_aux_vectors(u, v, w);
                    }
                }
            }
        }
        solver.addRows(lower_bounds.size(), &lower_bounds[0], &upper_bounds[0],  //
                       indices.size(), &starts[0], &indices[0], &values[0]);

        lower_bounds.clear();
        lower_bounds.shrink_to_fit();
        upper_bounds.clear();
        upper_bounds.shrink_to_fit();
        starts.clear();
        starts.shrink_to_fit();
        indices.clear();
        indices.shrink_to_fit();
        values.clear();
        values.shrink_to_fit();
    }

    /**
     * @brief sets all columns types to integral
     */
    inline void set_columns_integral() {
        std::vector<HighsVarType> integrality(get_n_cols(), HighsVarType::kInteger);
        solver.changeColsIntegrality(0, get_n_cols(), &integrality[0]);
    }
};

};  // namespace pace

#endif
