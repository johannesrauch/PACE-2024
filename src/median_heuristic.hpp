#ifndef PACE_MEDIAN_HEURISTIC_HPP
#define PACE_MEDIAN_HEURISTIC_HPP

#ifndef PACE_CONST_PROBMEDIAN_ITERATIONS
#define PACE_CONST_PROBMEDIAN_ITERATIONS 128
#endif

#include <vector>

#include "crossings.hpp"
#include "debug_printf.hpp"
#include "instance.hpp"
#include "random.hpp"
#include "shift_heuristic.hpp"
#include "vector_utils.hpp"

namespace pace {

/**
 * @brief median heuristic solver
 * based on Eades' and Wormald's paper https://doi.org/10.1007/BF01187020
 *
 * @tparam T vertex type
 * @tparam R crossing number type
 */
template <typename T, typename R>
class median_heuristic {
    const bipartite_graph<T>& graph;
    const std::size_t n_free;
    std::vector<T> medians;

   public:
    shift_heuristic<T, R> shift_h;

    /**
     * @brief initializes median heuristic for `instance`
     */
    median_heuristic(const instance<T, R>& instance)
        : graph(instance.graph()),  //
          n_free(graph.get_n_free()),
          medians(n_free),
          shift_h(instance) {
        fill_medians();
    }

    // delete copy and move constructor and copy and move assignment
    median_heuristic(const median_heuristic<T, R>& other) = delete;
    median_heuristic(median_heuristic<T, R>&& other) = delete;
    median_heuristic<T, R>& operator=(median_heuristic<T, R>& other) = delete;
    median_heuristic<T, R>& operator=(const median_heuristic<T, R>& other) = delete;

    /**
     * @brief compares the medians of a and b
     * we break ties by considering the parity of the degrees of a and b
     *
     * @return true a < b (according to medians)
     * @return false a >= b (according to medians)
     */
    bool compare(const T& a, const T& b) const {
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
     *
     * @param ordering out parameter where result is stored
     * @return uint32_t number of crossings
     */
    template <bool SHIFT = true>
    uint32_t operator()(std::vector<T>& ordering) {
        ordering.resize(n_free);
        for (T i = 0; i < n_free; ++i) ordering[i] = i;
        sort(ordering.begin(), ordering.end(),
             [=](const T& a, const T& b) -> bool { return this->compare(a, b); });
        if constexpr (SHIFT) {
            return shift_h(ordering);
        } else {
            return number_of_crossings(graph, ordering);
        }
    }

   private:
    /**
     * @brief computes medians of every adjacency_lists[i] using `median(...)`
     * and stores them in medians
     */
    inline void fill_medians() {
        for (std::size_t i = 0; i < n_free; ++i) {
            medians[i] = median(graph.get_neighbors_of_free(i));
        }
    }

    /**
     * @brief comparator using `medians`
     */
    struct compare {
        const std::size_t& n_free;
        const std::vector<T>& medians;

        /**
         * @brief compares the medians of a and b
         * we break ties by considering the parity of the degrees of a and b
         *
         * @return true a < b (according to medians)
         * @return false a >= b (according to medians)
         */
        bool operator()(const T& a, const T& b) const {
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
    };
};

/**
 * @brief probabilistic median heuristic solver
 * based on Nagamochi's paper https://doi.org/10.1007/s00454-005-1168-0
 *
 * @tparam T vertex type
 */
template <typename T, typename R>
class probmedian_heuristic {
    const bipartite_graph<T>& graph;
    const std::size_t n_free;
    std::vector<T> another_ordering;
    std::vector<T> medians;
    std::vector<T> randomized_medians;
    std::uniform_real_distribution<> distribution{0.0957, 0.9043};  // from Nagamochi's paper
    uint32_t lower_bound, upper_bound{0};
    median_heuristic<T, R> median_h;
    shift_heuristic<T, R> shift_h;

   public:
    /**
     * @brief construct a new probabilistic median heuristic object
     *
     * @param graph the instance
     * @param ordering vector where the ordering is stored
     * @param nof_iterations number of iterations
     */
    probmedian_heuristic(const instance<T, R>& instance)
        : graph(instance.graph()),
          n_free(graph.get_n_free()),
          another_ordering(n_free),
          medians(n_free),
          randomized_medians(n_free),
          lower_bound(instance.get_lower_bound()),
          median_h(instance),
          shift_h(instance) {
        for (T i = 0; i < n_free; ++i) another_ordering[i] = i;
        fill_medians();
    }

    // delete copy and move constructor and copy and move assignment
    probmedian_heuristic(const probmedian_heuristic<T, R>& other) = delete;
    probmedian_heuristic(probmedian_heuristic<T, R>&& other) = delete;
    probmedian_heuristic<T, R>& operator=(probmedian_heuristic<T, R>& other) = delete;
    probmedian_heuristic<T, R>& operator=(const probmedian_heuristic<T, R>& other) = delete;

    /**
     * @brief compares the randomized medians of a and b
     * we break ties by considering the medians of a and b
     * and then the parity of the degrees of a and b
     * (odd degrees are smaller)
     *
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

    /**
     * @brief runs the probabilistic median heuristic solver
     * and stores the result in ordering
     *
     * @return uint32_t number of crossings
     */
    template <std::size_t NOF_ITERATIONS = PACE_CONST_PROBMEDIAN_ITERATIONS>
    uint32_t operator()(std::vector<T>& ordering) {
        upper_bound = median_h(ordering);
        // try to find a better solution with probabilistic median heuristic
        for (std::size_t i = 0; lower_bound < upper_bound && i < NOF_ITERATIONS; ++i) {
            uint32_t candidate = generate_another_ordering();
            if (candidate < upper_bound) {
                upper_bound = candidate;
                ordering = another_ordering;
            }
        }
        return upper_bound;
    }

   private:
    inline uint32_t generate_another_ordering() {
        fill_randomized_medians();
        sort(another_ordering.begin(), another_ordering.end(),
             [=](const T& a, const T& b) -> bool { return this->compare(a, b); });
        return shift_h(another_ordering);
    }

    /**
     * @brief returns randomized median of `graph.get_neighbors_of_free(i)`
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
    inline void fill_medians() {
        for (std::size_t i = 0; i < n_free; ++i) {
            medians[i] = median(graph.get_neighbors_of_free(i));
        }
    }
};

};  // namespace pace

#endif