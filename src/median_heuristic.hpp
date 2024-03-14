#ifndef PACE2024_MEDIAN_HEURISTIC_HPP
#define PACE2024_MEDIAN_HEURISTIC_HPP

#include <vector>

#include "bipartite_graph.hpp"
#include "crossings.hpp"
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
    const bipartite_graph<T>& graph;
    const std::size_t n_free;
    std::vector<T>& ordering;
    std::vector<T> medians;

   public:
    /**
     * @brief construct a new median heuristic object
     *
     * @param graph the instance
     * @param ordering vector where we store the computed ordering
     */
    median_heuristic(const bipartite_graph<T>& graph, std::vector<T>& ordering)
        : graph(graph), n_free(graph.get_n_free()), ordering(ordering), medians(n_free) {
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
        assert(n_free == medians.size());
        assert(a < n_free && b < n_free);
        if (medians[a] < medians[b]) {
            return true;
        } else if (medians[a] > medians[b]) {
            return false;
        } else if (graph.degree_of_free(a) % 2 == 1 && graph.degree_of_free(b) % 2 == 0) {
            return true;
        } else if (graph.degree_of_free(a) % 2 == 0 && graph.degree_of_free(b) % 2 == 1) {
            return false;
        } else {
            return a < b;
        }
    }

    /**
     * @brief runs the median heuristic and stores the result in `ordering`
     */
    void run() {
        ordering.resize(n_free);
        for (T i = 0; i < n_free; ++i) ordering[i] = i;
        sort(ordering.begin(), ordering.end(),
             [=](const T& a, const T& b) -> bool { return this->compare(a, b); });
    }

   private:
    /**
     * @brief computes medians of every adjacency_lists[i] using `internal::median(...)`
     * and stores them in medians
     */
    void fill_medians() {
        for (std::size_t i = 0; i < n_free; ++i) {
            medians[i] = internal::median(graph.get_neighbors_of_free(i));
        }
    }
};

/**
 * @brief probabilistic median heuristic solver
 * based on Nagamochi's paper https://doi.org/10.1007/s00454-005-1168-0
 *
 * @tparam T vertex type
 */
template <typename T>
class probabilistic_median_heuristic {
   private:
    const bipartite_graph<T>& graph;
    const std::size_t n_free;
    std::vector<T>& ordering;
    std::vector<T> another_ordering;
    std::vector<T> medians;
    std::vector<T> randomized_medians;
    std::uniform_real_distribution<> distribution;
    uint32_t best;

   public:
    /**
     * @brief construct a new probabilistic median heuristic object
     *
     * @param graph the instance
     * @param ordering vector where the ordering is stored
     * @param nof_iterations number of iterations
     */
    probabilistic_median_heuristic(const bipartite_graph<T>& graph, std::vector<T>& ordering)
        : graph(graph),
          n_free(graph.get_n_free()),
          ordering(ordering),
          another_ordering(n_free),
          medians(n_free),
          randomized_medians(n_free),
          distribution(0.0957, 0.9043)  // from Nagamochi's paper
    {
        // compute initial solution with normal median heuristic
        median_heuristic<T>(graph, ordering).run();
        best = number_of_crossings(graph, ordering);

        // initialize vectors
        for (T i = 0; i < n_free; ++i) another_ordering[i] = i;
        fill_medians();
    }

    // delete copy and move constructor and copy and move assignment
    probabilistic_median_heuristic(const probabilistic_median_heuristic<T>& other) = delete;
    probabilistic_median_heuristic(probabilistic_median_heuristic<T>&& other) = delete;
    probabilistic_median_heuristic<T>& operator=(probabilistic_median_heuristic<T>& other) = delete;
    probabilistic_median_heuristic<T>& operator=(const probabilistic_median_heuristic<T>& other) =
        delete;

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
        assert(n_free == randomized_medians.size());
        assert(n_free == medians.size());
        assert(a < n_free && b < n_free);
        if (randomized_medians[a] < randomized_medians[b]) {
            return true;
        } else if (randomized_medians[a] > randomized_medians[b]) {
            return false;
        } else if (medians[a] < medians[b]) {  // from here random_medians[a] == random_medians[b]
            return true;
        } else if (medians[a] > medians[b]) {  // from here medians[a] == medians[b]
            return false;
        } else if (graph.degree_of_free(a) % 2 == 1 && graph.degree_of_free(b) % 2 == 0) {
            return true;
        } else if (graph.degree_of_free(a) % 2 == 0 && graph.degree_of_free(b) % 2 == 1) {
            return false;
        } else {
            return a < b;
        }
    }

    uint32_t get_best() { return best; }

    /**
     * @brief runs the probabilistic median heuristic solver
     * and stores the result in ordering
     *
     * @return R number of crossings
     */
    uint32_t run(const std::size_t nof_iterations = 10) {
        // try to find a better solution with probabilistic median heuristic
        for (std::size_t iteration = 0; iteration < nof_iterations; ++iteration) {
            // PACE2024_DEBUG_PRINTF("%s\n", iteration);
            fill_randomized_medians();
            sort(another_ordering.begin(), another_ordering.end(),
                 [=](const T& a, const T& b) -> bool { return this->compare(a, b); });

            uint32_t candidate = number_of_crossings(graph, another_ordering);
            if (candidate < best) {
                best = candidate;
                ordering = another_ordering;
            }
        }
        return best;
    }

   private:
    /**
     * @brief returns randomized median of `graph.get_neighbors_of_free(i)`
     *
     * @param i vertex
     * @return T randomized median
     */
    inline T randomized_median(const std::size_t& i) {
        const auto& neighbors = graph.get_neighbors_of_free(i);
        const std::size_t nof_neighbors = graph.degree_of_free(i);
        if (nof_neighbors == 0) {
            return 0;
        } else {
            const std::size_t j = static_cast<std::size_t>(distribution(generator) * nof_neighbors);
            assert(j < nof_neighbors);
            return neighbors[j];
        }
    }

    /**
     * @brief computes randomized medians of every `graph.get_neighbors_of_free(i)`
     * using `randomized_median(i)` and stores them in `randomized_medians`
     */
    void fill_randomized_medians() {
        for (std::size_t i = 0; i < n_free; ++i) {
            randomized_medians[i] = randomized_median(i);
        }
    }

    /**
     * @brief computes medians of every `graph.get_neighbors_of_free(i)` using
     * `internal::median(...)` and stores them in `medians`
     */
    void fill_medians() {
        for (std::size_t i = 0; i < n_free; ++i) {
            medians[i] = internal::median(graph.get_neighbors_of_free(i));
        }
    }
};

};  // namespace pace2024

#endif