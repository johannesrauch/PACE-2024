#ifndef PACE_HEURISTIC_RELAX_HEURISTIC_HPP
#define PACE_HEURISTIC_RELAX_HEURISTIC_HPP

#include "utils/crossings_utils.hpp"
#include "utils/restr_graph_utils.hpp"

namespace pace {

class relax_heuristic : public instance_view {
    uint8_t n_lookahead;
    std::size_t n_iterations{0};
    std::size_t n_cycles{0};
    std::vector<vertex_t> another_ordering;

   public:
    relax_heuristic(instance &instance_, std::uint8_t n_lookahead = 64)
        : instance_view(instance_), n_lookahead(n_lookahead) {}

    crossing_number_t operator()(const std::vector<double> &column_values, std::vector<vertex_t> &ordering) {
        build_restr_graph(column_values, unsettled_pairs(), restriction_graph());
        n_cycles = topological_sort_rd(restriction_graph(), ordering);
        crossing_number_t n_crossings = number_of_crossings(graph, ordering);

        n_iterations = 0;
        for (uint8_t i = 1; i < n_lookahead; ++i) {
            const crossing_number_t candidate = generate_another_ordering(column_values);
            if (candidate < n_crossings) {
                i = 0;
                n_crossings = candidate;
                std::swap(ordering, another_ordering);
            }
        }

        update_upper_bound(n_crossings);
        return n_crossings;
    }

   private:
    crossing_number_t generate_another_ordering(const std::vector<double> &column_values) {
        build_restr_graph(column_values, unsettled_pairs(), restriction_graph());
        n_cycles += topological_sort_rd(restriction_graph(), another_ordering);
        return number_of_crossings(graph, another_ordering);
    }
};

};  // namespace pace

#endif
