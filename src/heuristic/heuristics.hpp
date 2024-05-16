#ifndef PACE_HEURISTICS_HEURISTICS_HPP
#define PACE_HEURISTICS_HEURISTICS_HPP

#include "heuristic/barycenter_heuristic.hpp"
#include "heuristic/median_heuristic.hpp"
#include "heuristic/round_heuristic.hpp"
#include "heuristic/sort_heuristic.hpp"
#include "heuristic/rins_heuristic.hpp"

namespace pace {

class heuristics : public instance_view {
    std::vector<vertex_t> another_ordering;
    barycenter_heuristic barycenter_h;
    probmedian_heuristic probmedian_h;
    round_heuristic round_h;
    sort_heuristic sort_h;
    shift_heuristic shift_h;
    rins_heuristic rins_h;

   public:
    heuristics(instance &instance_)
        : instance_view(instance_),
          another_ordering(n_free),
          barycenter_h(instance_),
          probmedian_h(instance_),
          round_h(instance_),
          sort_h(instance_),
          shift_h(instance_),
          rins_h(instance_) {}

    crossing_number_t uninformed(std::vector<vertex_t> &ordering);

    crossing_number_t informed(highs_lp &lp, std::vector<vertex_t> &ordering, const bool do_rins);

    crossing_number_t shift(std::vector<vertex_t> &ordering,
                            const crossing_number_t n_cr);

    crossing_number_t shift(std::vector<vertex_t> &ordering) {
        return shift(ordering, number_of_crossings(graph, ordering));
    }

    double get_confidence() const { return round_h.get_confidence(); }
};

};  // namespace pace

#endif
