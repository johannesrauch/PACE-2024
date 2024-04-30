#ifndef PACE_HEURISTICS_HEURISTICS_HPP
#define PACE_HEURISTICS_HEURISTICS_HPP

#include "heuristic/barycenter_heuristic.hpp"
#include "heuristic/median_heuristic.hpp"
#include "heuristic/round_heuristic.hpp"
#include "heuristic/sort_heuristic.hpp"

namespace pace {

class heuristics : public instance_view {
    std::vector<vertex_t> another_ordering;
    barycenter_heuristic barycenter_h;
    probmedian_heuristic probmedian_h;
    round_heuristic round_h;
    sort_heuristic sort_h;

   public:
    heuristics(instance &instance_)
        : instance_view(instance_),
          another_ordering(n_free),
          barycenter_h(instance_),
          probmedian_h(instance_),
          round_h(instance_),
          sort_h(instance_) {}

    crossing_number_t uninformed(std::vector<vertex_t> &ordering) {
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

    crossing_number_t informed(highs_lp &lp, std::vector<vertex_t> &ordering) {
        PACE_DEBUG_PRINTF("start informed heuristics\n");
        crossing_number_t n_cr = sort_h(lp, ordering);
        if (n_cr <= lower_bound()) return n_cr;
        const crossing_number_t n_cr_r = round_h(lp, another_ordering);
        if (n_cr_r < n_cr) {
            n_cr = n_cr_r;
            std::swap(another_ordering, ordering);
        }
        PACE_DEBUG_PRINTF("end   informed heuristics, ub=%s\n", upper_bound);
        return n_cr;
    }

    double get_confidence() const { return round_h.get_confidence(); }
};

};  // namespace pace

#endif
