#ifndef PACE_EXACT_EXACT_SOLVER_HPP
#define PACE_EXACT_EXACT_SOLVER_HPP

#include "io/input.hpp"

namespace pace {

template <typename T, typename R>
class exact_solver {
    input<T, R> &in;

   public:
    exact_solver(input<T, R> &in) : in(in) { in.try_split(); }

    uint32_t operator()(std::vector<T> &ordering) {
        if (in.get_n_subgraphs() >= 2) {
        }
    }

   private:
    uint32_t solve_if_split(std::vector<T> &ordering_) {
        const std::size_t n_subgraphs = in.get_n_subgraphs();
        assert(n_subgraphs >= 2);
        const bool first_graph_empty = in.is_first_graph_empty();
        std::vector<std::vector<T>> orderings;

        if (first_graph_empty) {
            const std::size_t n_free = in.get_subgraph(0).get_n_free();
            orderings.emplace_back(n_free);
            std::vector<T> &ordering = *orderings.rbegin();
            for (std::size_t i = 0; i < n_free) ordering[i] = i;
        }

        for (std::size_t i = first_graph_empty; i < n_subgraphs; ++i) {
            const std::size_t n_free = in.get_subgraph(0).get_n_free();
            orderings.emplace_back(n_free);
        }
    }
};

};  // namespace pace

#endif
