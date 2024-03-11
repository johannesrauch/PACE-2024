#ifndef PACE2024_MEDIAN_HEURISTIC_HPP
#define PACE2024_MEDIAN_HEURISTIC_HPP

#include <vector>

#include "bipartite_graph.hpp"
#include "crossing_number.hpp"
#include "debug_printf.hpp"
#include "random.hpp"

namespace pace2024 {

namespace internal {

/**
 * @brief returns the median of the vector vec
 *
 * @tparam T element type of vec
 * @param vec vector
 * @return T the median of vec
 */
template <typename T>
inline T median(const std::vector<T>& vec) {
    const std::size_t len = vec.size();
    if (len == 0) {
        return 0;
    } else if (len == 1) {
        return vec[0];
    } else {
        const std::size_t len2 = len / 2;
        return len % 2 == 0 ? (vec[len2 - 1] + vec[len2]) / 2 : vec[len / 2];
    }
}

};  // namespace internal

/**
 * @brief median heuristic solver
 * based on Eades' and Wormald's paper https://doi.org/10.1007/BF01187020
 *
 * @tparam T vertex type of the instance
 */
template <typename T>
class median_heuristic {
   private:
    const bipartite_graph<T>& input_graph;
    const std::vector<std::vector<T>>& adjacency_lists;
    const std::size_t n1;
    std::vector<T>& ordering;
    std::vector<T> medians;

   public:
    /**
     * @brief construct a new median heuristic object
     *
     * @param graph the instance
     * @param ordering vector where we store the computed ordering
     */
    median_heuristic(const bipartite_graph<T>& input_graph,
                     std::vector<T>& ordering)
        : input_graph(input_graph),
          adjacency_lists(input_graph.get_adjacency_lists()),
          n1(input_graph.get_n_free()),
          ordering(ordering),
          medians(n1) {
        fill_medians();
    }

    // delete copy and move constructor and copy and move assignment
    median_heuristic(const median_heuristic<T>& other) = delete;
    median_heuristic(median_heuristic<T>&& other) = delete;
    median_heuristic<T>& operator=(median_heuristic<T>& other) = delete;
    median_heuristic<T>& operator=(const median_heuristic<T>& other) = delete;

    /**
     * @brief compares the medians of a and b
     * we break ties by considering the parity of the degrees of a and b
     * if both degrees are odd, we flip a coin
     *
     * @param a
     * @param b
     * @return true a < b (according to medians)
     * @return false a >= b (according to medians)
     */
    bool compare(const T& a, const T& b) const {
        assert(n1 == medians.size());
        assert(a < n1 && b < n1);
        if (medians[a] < medians[b]) {
            return true;
        } else if (medians[a] > medians[b]) {
            return false;
        } else if (adjacency_lists[a].size() % 2 == 1 && adjacency_lists[b].size() % 2 == 0) {
            return true;
        } else if (adjacency_lists[a].size() % 2 == 0 && adjacency_lists[b].size() % 2 == 1) {
            return false;
        } else {
            return a < b;
        }
    }

    /**
     * @brief runs the median heuristic and stores the result in ordering
     */
    void run() {
        ordering.resize(n1);
        for (T i = 0; i < n1; ++i) ordering[i] = i;
        sort(ordering.begin(), ordering.end(), [=](const T& a, const T& b) -> bool {
            return this->compare(a, b);
        });
    }

   private:
    /**
     * @brief computes medians of every adjacency_lists[i] using internal::compute_median
     * and stores them in medians
     */
    void fill_medians() {
        for (std::size_t i = 0; i < n1; ++i) {
            medians[i] = internal::median(adjacency_lists[i]);
        }
    }
};

/**
 * @brief probabilistic median heuristic solver
 * based on Nagamochi's paper https://doi.org/10.1007/s00454-005-1168-0
 *
 * @tparam T
 * @tparam R
 */
template <typename T, typename R>
class probabilistic_median_heuristic {
   private:
    const bipartite_graph<T>& input_graph;
    const std::vector<std::vector<T>>& adjacency_lists;
    const std::size_t n1;
    std::vector<T>& ordering;
    std::vector<T> another_ordering;
    std::vector<T> medians;
    std::vector<T> randomized_medians;
    std::uniform_real_distribution<> distribution;
    R best;

   public:
    /**
     * @brief construct a new probabilistic median heuristic object
     *
     * @param input_graph the instance
     * @param ordering vector where the ordering is stored
     * @param nof_iterations number of iterations
     */
    probabilistic_median_heuristic(const bipartite_graph<T>& input_graph,
                                   std::vector<T>& ordering)
        : input_graph(input_graph),
          adjacency_lists(input_graph.get_adjacency_lists()),
          n1(input_graph.get_n_free()),
          ordering(ordering),
          another_ordering(n1),
          medians(n1),
          randomized_medians(n1),
          distribution(0.0957, 0.9043)  // from Nagamochi's paper
    {
        // compute initial solution with normal median heuristic
        median_heuristic<T>(input_graph, ordering).run();
        best = number_of_crossings<T, R>(input_graph, ordering);

        // initialize vectors
        for (T i = 0; i < n1; ++i) another_ordering[i] = i;
        fill_medians();
    }

    // delete copy and move constructor and copy and move assignment
    probabilistic_median_heuristic(const probabilistic_median_heuristic<T, R>& other) = delete;
    probabilistic_median_heuristic(probabilistic_median_heuristic<T, R>&& other) = delete;
    probabilistic_median_heuristic<T, R>& operator=(probabilistic_median_heuristic<T, R>& other) = delete;
    probabilistic_median_heuristic<T, R>& operator=(const probabilistic_median_heuristic<T, R>& other) = delete;

    /**
     * @brief compares the randomized medians of a and b
     * we break ties by considering the medians of a and b
     * and then the parity of the degrees of a and b
     * (odd degrees are smaller)
     *
     * @param a
     * @param b
     * @return true a < b (according to medians)
     * @return false a >= b (according to medians)
     */
    bool compare(const T& a, const T& b) const {
        assert(n1 == randomized_medians.size());
        assert(n1 == medians.size());
        assert(a < n1 && b < n1);
        if (randomized_medians[a] < randomized_medians[b]) {
            return true;
        } else if (randomized_medians[a] > randomized_medians[b]) {
            return false;
        } else if (medians[a] < medians[b]) {  // from here random_medians[a] == random_medians[b]
            return true;
        } else if (medians[a] > medians[b]) {  // from here medians[a] == medians[b]
            return false;
        } else if (adjacency_lists[a].size() % 2 == 1 && adjacency_lists[b].size() % 2 == 0) {
            return true;
        } else if (adjacency_lists[a].size() % 2 == 0 && adjacency_lists[b].size() % 2 == 1) {
            return false;
        } else {
            return a < b;
        }
    }

    R get_best() {
        return best;
    }

    /**
     * @brief runs the probabilistic median heuristic solver
     * and stores the result in ordering
     *
     * @return R number of crossings
     */
    R run(const std::size_t nof_iterations = 10) {
        // try to find a better solution with probabilistic median heuristic
        for (std::size_t iteration = 0; iteration < nof_iterations; ++iteration) {
            // PACE2024_DEBUG_PRINTF("%s\n", iteration);
            fill_randomized_medians();
            sort(another_ordering.begin(), another_ordering.end(), [=](const T& a, const T& b) -> bool {
                return this->compare(a, b);
            });

            R candidate = number_of_crossings<T, R>(input_graph, another_ordering);
            if (candidate < best) {
                best = candidate;
                ordering = another_ordering;
            }
        }
        return best;
    }

   private:
    /**
     * @brief returns randomized median of adjacency_lists[i]
     *
     * @param i vertex
     * @return T randomized median
     */
    inline T randomized_median(const std::size_t& i) {
        std::size_t nof_neighbors = adjacency_lists[i].size();
        if (nof_neighbors == 0) {
            return 0;
        } else {
            const std::size_t j = static_cast<std::size_t>(distribution(generator) * nof_neighbors);
            assert(j < nof_neighbors);
            return adjacency_lists[i][j];
        }
    }

    /**
     * @brief computes randomized medians of every adjacency_lists[i] using compute_randomized_median
     * and stores them in randomized_medians
     */
    void fill_randomized_medians() {
        for (std::size_t i = 0; i < n1; ++i) {
            randomized_medians[i] = randomized_median(i);
        }
    }

    /**
     * @brief computes medians of every adjacency_lists[i] using internal::compute_median
     * and stores them in medians
     */
    void fill_medians() {
        for (std::size_t i = 0; i < n1; ++i) {
            medians[i] = internal::median(adjacency_lists[i]);
        }
    }
};

};  // namespace pace2024

#endif