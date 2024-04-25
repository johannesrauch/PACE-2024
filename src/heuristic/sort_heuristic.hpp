#ifndef PACE_HEURISTIC_SORT_HEURISTIC_HPP
#define PACE_HEURISTIC_SORT_HEURISTIC_HPP

#include "exact/highs_lp.hpp"
#include "heuristic/shift_heuristic.hpp"
#include "utils/vector_utils.hpp"

namespace pace {

class sort_heuristic : public instance_view {
    highs_lp &lp;
    shift_heuristic shift_h;

   public:
    sort_heuristic(instance &instance_, highs_lp &lp) : instance_view(instance_), lp(lp), shift_h(instance_) {}

    crossing_number_t operator()(std::vector<vertex_t> &ordering) {
        std::vector<double> precedence(n_free);
        for (vertex_t u = 0; u < n_free; ++u) {
            for (vertex_t v = 0; v < u; ++v)
                precedence[u] += lp_solver_ptr->get_variable_value(flat_index(n_free, n_free_2, v, u));
            for (vertex_t v = u + 1u; v < n_free; ++v)
                precedence[u] += 1. - lp_solver_ptr->get_variable_value(flat_index(n_free, n_free_2, u, v));
        }
        identity(n_free, ordering);
        std::sort(ordering.begin(), ordering.end(),
                  [&](const vertex_t &u, const vertex_t &v) -> bool { return precedence[u] < precedence[v]; });
        return shift_h(ordering);
    }
};

};  // namespace pace

#endif
