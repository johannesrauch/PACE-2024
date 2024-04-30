#ifndef PACE_HEURISTIC_RELAX_HEURISTIC_HPP
#define PACE_HEURISTIC_RELAX_HEURISTIC_HPP

#include "heuristic/shift_heuristic.hpp"
#include "utils/crossings_utils.hpp"
#include "utils/restr_graph_utils.hpp"

namespace pace {

struct relax_heuristic_params {
    const bool do_shift{true};
    const uint8_t n_lookahead{16};
};

class relax_heuristic : public instance_view {
    relax_heuristic_params params;
    std::size_t n_restr_graphs_generated{0};
    std::size_t n_cycles{0};
    std::vector<vertex_t> another_ordering;
    shift_heuristic shift_h;

   public:
    relax_heuristic(instance &instance_, relax_heuristic_params params = relax_heuristic_params())
        : instance_view(instance_), params(params), shift_h(instance_) {}

    crossing_number_t operator()(const std::vector<double> &column_values, std::vector<vertex_t> &ordering) {
        n_restr_graphs_generated = 0;
        n_cycles = 0;
        crossing_number_t n_crossings = generate_ordering(column_values, ordering);

        for (uint8_t i = 1; i < params.n_lookahead; ++i) {
            const crossing_number_t candidate = generate_ordering(column_values, another_ordering);
            if (candidate < n_crossings) {
                i = 0;
                n_crossings = candidate;
                std::swap(ordering, another_ordering);
            }
        }

        update_ordering(ordering, n_crossings);
        return n_crossings;
    }

    double get_confidence() const {
        if (n_restr_graphs_generated == 0) return -1.;
        return static_cast<double>(n_cycles) / n_restr_graphs_generated;
    }

   private:
    crossing_number_t  //
    generate_ordering(const std::vector<double> &column_values, std::vector<vertex_t> &ordering) {
        build_restr_graph(column_values, unsettled_pairs(), restriction_graph());
        n_cycles += topological_sort_rd(restriction_graph(), ordering);
        ++n_restr_graphs_generated;
        if (params.do_shift)
            return shift_h(ordering);
        else
            return number_of_crossings(graph, ordering);
    }
};

};  // namespace pace

#endif
