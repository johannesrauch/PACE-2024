#include "heuristic/heuristics.hpp"

namespace pace {

crossing_number_t heuristics::uninformed(std::vector<vertex_t> &ordering) {
    PACE_DEBUG_PRINTF("start uninformed heuristics\n");
    crossing_number_t n_cr = barycenter_h(ordering);
    if (n_cr <= lower_bound()) return n_cr;
    const crossing_number_t n_cr_p = probmedian_h(another_ordering);
    if (n_cr_p < n_cr) {
        n_cr = n_cr_p;
        std::swap(another_ordering, ordering);
    }
    PACE_DEBUG_PRINTF("end   uninformed heuristics\n");
    return n_cr;
}

crossing_number_t heuristics::informed(  //
    highs_lp &lp, std::vector<vertex_t> &ordering, const bool do_rins) {
    PACE_DEBUG_PRINTF("start informed heuristics\n");

    crossing_number_t n_cr = sort_h(lp, ordering);
    const crossing_number_t n_cr_round = round_h(lp, another_ordering);
    if (n_cr_round < n_cr) {
        n_cr = n_cr_round;
        std::swap(another_ordering, ordering);
    }
    if (do_rins) {
        const crossing_number_t n_cr_rins = rins_h(lp, another_ordering);
        if (n_cr_rins < n_cr) {
            n_cr = n_cr_rins;
            std::swap(another_ordering, ordering);
        }
    }

    PACE_DEBUG_PRINTF("end   informed heuristics, ub=%s\n", upper_bound);
    return n_cr;
}

crossing_number_t heuristics::shift(  //
    std::vector<vertex_t> &ordering, const crossing_number_t n_cr) {
    return shift_h(ordering, n_cr);
}

}  // namespace pace
