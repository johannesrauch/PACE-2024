#include "heuristic/sort_heuristic.hpp"

#include "exact/highs_lp.hpp"

namespace pace {

crossing_number_t sort_heuristic::operator()(  //
    highs_lp &lp, std::vector<vertex_t> &ordering) {
    std::vector<double> precedence(n_free);
    for (vertex_t u = 0; u < n_free; ++u) {
        for (vertex_t v = 0; v < u; ++v)
            precedence[u] +=
                lp.get_variable_value(flat_index(n_free, n_free_2, v, u));
        for (vertex_t v = u + 1u; v < n_free; ++v)
            precedence[u] +=
                1. - lp.get_variable_value(flat_index(n_free, n_free_2, u, v));
    }
    identity(n_free, ordering);
    std::sort(ordering.begin(), ordering.end(),
              [&](const vertex_t &u, const vertex_t &v) -> bool {
                  return precedence[u] < precedence[v];
              });
    return shift_h(ordering);
}

}  // namespace pace