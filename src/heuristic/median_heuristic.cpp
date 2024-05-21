#include "heuristic/median_heuristic.hpp"

namespace pace {

crossing_number_t median_heuristic::operator()(
    std::vector<vertex_t>& ordering) {
    identity(n_free, ordering);
    sort(ordering.begin(), ordering.end(),
         [=](const vertex_t& a, const vertex_t& b) -> bool {
             return this->compare(a, b);
         });
    return shift_h(ordering);
}

crossing_number_t probmedian_heuristic::operator()(
    std::vector<vertex_t>& ordering, crossing_number_t n_crossings) {
    PACE_DEBUG_PRINTF("start probmedian heuristic\n");
    // try to find a better solution with probabilistic median heuristic
    for (uint8_t i = 1;
         lower_bound() < n_crossings && i <= params.n_lookahead &&
         !timelimit::was_sigterm_sent();
         ++i) {
        const crossing_number_t candidate = generate_another_ordering();
        if (candidate < n_crossings) {
            i = 0;
            n_crossings = candidate;
            std::swap(ordering, another_ordering);
        }
        PACE_DEBUG_PRINTF("%11s=%11u, %11s=%11u\n", "i", i, "pm look",
                          params.n_lookahead);
    }
    update_ordering(ordering, n_crossings);
    PACE_DEBUG_PRINTF("end   probmedian heuristic\n");
    return n_crossings;
}

crossing_number_t probmedian_heuristic::operator()(
    std::vector<vertex_t>& ordering) {
    const crossing_number_t n_crossings = median_h(ordering);
    // try to find a better solution with probabilistic median heuristic
    return (*this)(ordering, n_crossings);
}

}  // namespace pace