#ifndef PACE2024_MEDIAN_HEURISTIC_HPP
#define PACE2024_MEDIAN_HEURISTIC_HPP

#include <random>
#include <vector>

#include "bipartite_graph.hpp"
#include "crossings.hpp"

namespace pace2024 {

/**
 * @brief computes a linear ordering with the median heuristic
 *
 * @tparam T an unsigned integer type
 * @param graph the bipartite graph instance
 * @param ordering an empty array serving as output
 */
template <typename T>
void median_heuristic(const general_bipartite_graph<T>& graph, std::vector<T>& ordering) {
    const T n1 = graph.get_n1();
    ordering.resize(n1);
    std::vector<T> medians(n1);
    auto adjacency_lists = graph.get_adjacency_lists();

    for (std::size_t i = 0; i < n1; ++i) {
        ordering[i] = i;

        std::size_t nof_neighbors = adjacency_lists[i].size();
        if (nof_neighbors == 0) {
            medians[i] = 0;
        } else {
            medians[i] = adjacency_lists[i][nof_neighbors / 2];
        }
    }

    sort(ordering.begin(), ordering.end(), [&](const T& a, const T& b) -> bool {
        return medians[a] < medians[b];
    });
}

std::mt19937 generator(std::random_device{}());
std::uniform_real_distribution<> distribution(0.0957, 0.9043);

/**
 * @brief
 *
 * @tparam T
 * @param graph
 * @param ordering
 */
template <typename T>
void prob_median_heuristic(const general_bipartite_graph<T>& graph,
                           const folded_square_matrix<T> cr_matrix,
                           std::vector<T>& ordering,
                           std::size_t nof_iterations = 1000) {
    // compute a solution with the normal median heuristic
    median_heuristic(graph, ordering);
    T best = compute_crossings(graph, cr_matrix, ordering);

    const T n1 = graph.get_n1();
    std::vector<T> another_ordering(n1);
    for (std::size_t i = 0; i < n1; ++i) {
            another_ordering[i] = i;
    }
    std::vector<T> medians(n1);
    auto adjacency_lists = graph.get_adjacency_lists();

    // run the probabilistic algorithm to find a better solution
    for (std::size_t iteration = 0; iteration < nof_iterations; ++iteration) {
        for (std::size_t i = 0; i < n1; ++i) {
            std::size_t nof_neighbors = adjacency_lists[i].size();
            if (nof_neighbors == 0) {
                medians[i] = 0;
            } else {
                medians[i] = adjacency_lists[i][(std::size_t)(distribution(generator) * nof_neighbors)];
            }
        }

        sort(another_ordering.begin(), another_ordering.end(), [&](const T& a, const T& b) -> bool {
            return medians[a] < medians[b];
        });

        T candidate = compute_crossings(graph, cr_matrix, another_ordering);
        if (candidate < best) {
            best = candidate;
            ordering = another_ordering;
        }
    }
}

};  // namespace pace2024

#endif