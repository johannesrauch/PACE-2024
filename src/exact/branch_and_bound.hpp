#ifndef PACE_EXACT_BRANCH_AND_BOUND_HPP
#define PACE_EXACT_BRANCH_AND_BOUND_HPP

#include <map>

#include "crossings.hpp"
#include "digraph.hpp"
#include "heuristics.hpp"
#include "instance.hpp"
#include "oracle.hpp"
#include "topological_sort.hpp"

namespace pace {

template <typename T, typename R>
class branch_and_bound {
    const pace::instance<T, R> &instance;
    const bipartite_graph<T> &graph;
    const folded_matrix<R> &cr_matrix;
    const std::size_t n_free;

    std::vector<std::pair<T, T>> unsettled;
    std::vector<bool> settled;
    digraph<T> restriction_graph;

    std::vector<T> ordering;
    const uint32_t &lower_bound;
    uint32_t upper_bound;

   public:
    branch_and_bound(const pace::instance<T, R> &instance)
        : instance(instance),
          graph(instance.graph()),
          cr_matrix(instance.cr_matrix()),
          n_free(graph.get_n_free()),
          restriction_graph(n_free),
          lower_bound(instance.get_lower_bound()) {}

    bool compare(const std::pair<T, T> &p, const std::pair<T, T> &q) {
        const uint32_t c_p = cr_matrix(p.first, p.second) + cr_matrix(p.second, p.first);
        const uint32_t c_q = cr_matrix(q.first, q.second) + cr_matrix(q.second, q.first);
        if (c_p > c_q)
            return true;
        else if (c_p < c_q)
            return false;
        else
            return p < q;
    }

    uint32_t operator()() {
        PACE_DEBUG_PRINTF("\n\nstart heuristic\n");
        upper_bound = heuristics(instance, ordering);
        PACE_DEBUG_PRINTF("end   heuristic\n");

        uint32_t objective_offset = 0;
        {
            PACE_DEBUG_PRINTF("start oracle\n");
            std::vector<int> magic;
            oracle<T, R> oracle(instance, upper_bound);
            objective_offset = oracle.build(restriction_graph, magic, unsettled);
            PACE_DEBUG_PRINTF("end   oracle\n");
        }
        settled.clear();
        settled.resize(unsettled.size());

        std::sort(unsettled.begin(), unsettled.end(),
                  [=](const std::pair<T, T> &p, const std::pair<T, T> &q) -> bool { return this->compare(p, q); });

        branch(0, objective_offset);
        return upper_bound;
    }

    const std::vector<T> &get_ordering() { return ordering; }

   private:
    void branch(std::size_t j, const uint32_t objective_value) {
        if (objective_value >= upper_bound) return;
        while (j < settled.size() && settled[j]) ++j;
        if (j >= unsettled.size()) {
            upper_bound = objective_value;
            topological_sort(restriction_graph, ordering);
            assert(upper_bound == number_of_crossings(graph, ordering));
            return;
        }

        const auto &[u, v] = unsettled[j];
        std::map<T, T> old_degrees;
        std::vector<std::pair<T, T>> new_arcs;

        uint32_t objective_gain = fix(u, v, old_degrees, new_arcs, j);
        branch(j + 1, objective_gain + objective_value);
        unfix(old_degrees, new_arcs, j);

        objective_gain = fix(v, u, old_degrees, new_arcs, j);
        branch(j + 1, objective_gain + objective_value);
        unfix(old_degrees, new_arcs, j);
    }

    std::size_t find(const T &u, const T &v) {
        const std::pair<T, T> p = u < v ? std::make_pair(u, v) : std::make_pair(v, u);
        auto it = std::lower_bound(
            unsettled.begin(), unsettled.end(), p,
            [=](const std::pair<T, T> &p, const std::pair<T, T> &q) -> bool { return this->compare(p, q); });
        assert(*it == p);
        return std::distance(unsettled.begin(), it);
    }

    uint32_t fix(const T &u, const T &v, std::map<T, T> &old_degrees, std::vector<std::pair<T, T>> &new_arcs,
                 const std::size_t &j) {
        old_degrees.clear();
        new_arcs.clear();
        uint32_t objective_gain = 0;

        old_degrees[u] = restriction_graph.degree(u);
        restriction_graph.add_arc(u, v);
        objective_gain += cr_matrix(u, v);
        settled[j] = true;

        transitive_hull(restriction_graph, new_arcs);
        for (const auto &[x, y] : new_arcs) {
            if (old_degrees.count(x) == 0) {
                old_degrees[x] = restriction_graph.degree(x);
            }
            restriction_graph.add_arc(x, y);
            objective_gain += cr_matrix(x, y);
            settled[find(x, y)] = true;
        }

        return objective_gain;
    }

    void unfix(const std::map<T, T> &old_degrees, const std::vector<std::pair<T, T>> &new_arcs, const std::size_t &j) {
        settled[j] = false;
        for (const auto &[x, d_x] : old_degrees) {
            restriction_graph.resize_neighbors(x, d_x);
        }
        for (const auto &[x, y] : new_arcs) {
            settled[find(x, y)] = false;
        }
    }
};

};  // namespace pace

#endif
