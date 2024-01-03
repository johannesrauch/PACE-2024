#ifndef PACE2024_MEDIAN_HEURISTIC_HPP
#define PACE2024_MEDIAN_HEURISTIC_HPP

#include <random>
#include <vector>

#include "bipartite_graph.hpp"
#include "crossings.hpp"
#include "debug.hpp"
#include "random.hpp"

namespace pace2024 {

/**
 * @brief sets ordering[i] = i for every i
 *
 * @tparam T
 * @param ordering
 */
template <typename T>
inline void fill_with_identity(std::vector<T>& ordering) {
    std::size_t n1 = ordering.size();
    for (T i = 0; i < n1; ++i) {
        ordering[i] = i;
    }
}

/**
 * @brief computes and returns the median of adjacency_list
 * if size = adjacency_list.size() is odd, it is adjacency_list[size / 2]
 * otherwise it is (adjacency_list[size / 2 - 1] + adjacency_list[size / 2]) / 2
 *
 * @tparam T
 * @param adjacency_list
 * @return T median of adjacency_list
 */
template <typename T>
inline T compute_median(const std::vector<T>& adjacency_list) {
    const std::size_t nof_neighbors = adjacency_list.size();
    if (nof_neighbors == 0) {
        return 0;
    } else {
        // todo: median for nof_neighbors even
        return adjacency_list[nof_neighbors / 2];
    }
}

/**
 * @brief computes medians of every adjacency_lists[i] using compute_median
 * and stores them in medians
 *
 * @tparam T
 * @param adjacency_lists
 * @param medians
 */
template <typename T>
inline void compute_medians(const general_bipartite_graph<T>& graph, std::vector<T>& medians) {
    const std::size_t n1 = graph.get_n1();
    auto adjacency_lists = graph.get_adjacency_lists();
    medians.resize(n1);
    for (std::size_t i = 0; i < n1; ++i) {
        const T median = compute_median(adjacency_lists[i]);
        assert(median < graph.get_n0());
        medians[i] = median;
    }
}

template <typename T>
struct median_heuristic_comparator {
    const std::vector<std::vector<T>>& adjacency_lists;
    const std::vector<T>& medians;

    median_heuristic_comparator(const std::vector<std::vector<T>>& adjacency_lists,
                                const std::vector<T>& medians)
        : adjacency_lists(adjacency_lists),
          medians(medians) {
        assert(adjacency_lists.size() == medians.size());
    }

    bool operator()(const T& a, const T& b) {
        assert(a < adjacency_lists.size() && b < adjacency_lists.size());
        if (medians[a] < medians[b]) {
            return true;
        } else if (medians[a] > medians[b]) {
            return false;
        } else if (adjacency_lists[a].size() % 2 == 1) {  // from here medians[a] == medians[b]
            return true;
        } else if (adjacency_lists[b].size() % 2 == 1) {
            return false;
        } else {
            return coinflip();
        }
    }
};

/**
 * @brief computes a linear ordering with the median heuristic
 *
 * @tparam T an unsigned integer type
 * @param graph the bipartite input graph
 * @param ordering an empty array serving as output
 */
template <typename T>
void median_heuristic(const general_bipartite_graph<T>& graph,
                      std::vector<T>& ordering) {
    const std::size_t n1 = graph.get_n1();

    ordering.resize(n1);
    fill_with_identity(ordering);

    std::vector<T> medians(n1);
    const std::vector<std::vector<T>>& adjacency_lists = graph.get_adjacency_lists();
    compute_medians(graph, medians);

    median_heuristic_comparator<T> comparator(adjacency_lists, medians);
    sort(ordering.begin(), ordering.end(), comparator);
}

/**
 * @brief computes and returns a randomized "median" of adjacency_list
 * draws a pseudorandom value r from the distribution
 * returns adjacency_list[floor(r * adjacency_list.size())]
 *
 * @tparam T
 * @param adjacency_list
 * @param distribution a distribution producing values in [0, 1]
 * @return T a randomized median of adjacency_list
 */
template <typename T>
inline T compute_random_median(const std::vector<T>& adjacency_list, std::uniform_real_distribution<>& distribution) {
    std::size_t nof_neighbors = adjacency_list.size();
    if (nof_neighbors == 0) {
        return 0;
    } else {
        return adjacency_list[(std::size_t)(distribution(generator) * nof_neighbors)];
    }
}

template <typename T>
struct prob_median_heuristic_comparator {
    const std::vector<std::vector<T>>& adjacency_lists;
    const std::vector<T>& random_medians;
    const std::vector<T>& medians;

    prob_median_heuristic_comparator(const std::vector<std::vector<T>>& adjacency_lists,
                                     const std::vector<T>& random_medians,
                                     const std::vector<T>& medians)
        : adjacency_lists(adjacency_lists),
          random_medians(random_medians),
          medians(medians) {
        assert(adjacency_lists.size() == random_medians.size() && random_medians.size() == medians.size());
    }

    bool operator()(const T& a, const T& b) {
        assert(a < medians.size() && b < medians.size());
        if (random_medians[a] < random_medians[b]) {
            return true;
        } else if (random_medians[a] > random_medians[b]) {
            return false;
        } else if (medians[a] < medians[b]) {  // from here random_medians[a] == random_medians[b]
            return true;
        } else if (medians[a] > medians[b]) {  // from here medians[a] == medians[b]
            return false;
        } else if (adjacency_lists[a].size() % 2 == 1) {
            return true;
        } else if (adjacency_lists[b].size() % 2 == 1) {
            return false;
        } else {
            return coinflip();
        }
    }
};

/**
 * @brief first calls median_heuristic for an initial solution
 * after that, we try to find a better solution with a randomized median
 * heuristic
 *
 * @tparam T an unsigned integer type
 * @param graph
 * @param ordering an empty array serving as output
 */
template <typename T>
T prob_median_heuristic(const general_bipartite_graph<T>& graph,
                        const folded_square_matrix<T>& cr_matrix,
                        std::vector<T>& ordering,
                        std::size_t nof_iterations = 1000) {
    // compute a solution with the normal median heuristic
    median_heuristic(graph, ordering);
    T best = compute_crossings(cr_matrix, ordering);

    const T n1 = graph.get_n1();
    std::vector<T> another_ordering(n1);
    // another_ordering[i] = i for all i
    fill_with_identity(another_ordering);

    std::vector<T> medians(n1);
    compute_medians(medians);

    std::uniform_real_distribution<> distribution(0.0957, 0.9043);
    std::vector<T> random_medians(n1);
    auto adjacency_lists = graph.get_adjacency_lists();

    // run the probabilistic algorithm to find a better solution
    prob_median_heuristic_comparator<T> comparator(adjacency_lists, random_medians, medians);
    for (std::size_t iteration = 0; iteration < nof_iterations; ++iteration) {
        for (std::size_t i = 0; i < n1; ++i) {
            random_medians[i] = compute_random_median(adjacency_lists[i], distribution);
        }

        sort(another_ordering.begin(), another_ordering.end(), comparator);

        T candidate = compute_crossings(cr_matrix, another_ordering);
        if (candidate < best) {
            best = candidate;
            ordering = another_ordering;
        }
    }

    return best;
}

};  // namespace pace2024

#endif