#ifndef PACE_HEURISTICS_MEDIAN_HEURISTIC_HPP
#define PACE_HEURISTICS_MEDIAN_HEURISTIC_HPP

#include "heuristic/shift_heuristic.hpp"
#include "utils/randomness_utils.hpp"
#include "utils/vector_utils.hpp"

namespace pace {

struct median_heuristic_params {
    const uint8_t n_lookahead{
#ifndef NDEBUG
        16
#else
        64
#endif
    };
};

/**
 * @brief median heuristic solver
 * based on Eades' and Wormald's paper https://doi.org/10.1007/BF01187020
 */
class median_heuristic : public instance_view {
    std::vector<vertex_t> medians;
    shift_heuristic shift_h;

   public:
    /**
     * @brief initializes median heuristic for `instance`
     */
    median_heuristic(instance& instance_)
        : instance_view(instance_), medians(n_free), shift_h(instance_) {
        fill_medians();
    }

    // delete copy and move constructor and copy and move assignment
    median_heuristic(const median_heuristic& other) = delete;
    median_heuristic(median_heuristic&& other) = delete;
    median_heuristic& operator=(median_heuristic& other) = delete;
    median_heuristic& operator=(const median_heuristic& other) = delete;

    /**
     * @brief runs the median heuristic and stores the result in 'ordering'
     *
     * @param ordering out parameter where result is stored
     * @return crossing_number_t number of crossings of 'ordering'
     */
    crossing_number_t operator()(std::vector<vertex_t>& ordering) {
        identity(n_free, ordering);
        sort(ordering.begin(), ordering.end(),
             [=](const vertex_t& a, const vertex_t& b) -> bool {
                 return this->compare(a, b);
             });
        return shift_h(ordering);
    }

    shift_heuristic& get_shift_heuristic() { return shift_h; }

    /**
     * @brief compares the medians of a and b
     * we break ties by considering the parity of the degrees of a and b
     *
     * @return true a < b (according to medians)
     * @return false a >= b (according to medians)
     */
    inline bool compare(const vertex_t& a, const vertex_t& b) const {
        assert(a < n_free && b < n_free);
        if (medians[a] < medians[b])
            return true;
        else if (medians[a] > medians[b])
            return false;
        else if (graph.get_degree(a) % 2 == 1 && graph.get_degree(b) % 2 == 0)
            return true;
        else if (graph.get_degree(a) % 2 == 0 && graph.get_degree(b) % 2 == 1)
            return false;
        else
            return a < b;
    }

   private:
    /**
     * @brief computes medians of every adjacency_lists[i] using `median(...)`
     * and stores them in medians
     */
    inline void fill_medians() {
        for (std::size_t i = 0; i < n_free; ++i) {
            medians[i] = median(graph.get_neighbors(i));
        }
    }
};

/**
 * @brief probabilistic median heuristic solver
 * based on Nagamochi's paper https://doi.org/10.1007/s00454-005-1168-0
 */
class probmedian_heuristic : public instance_view {
    std::vector<vertex_t> another_ordering;
    std::vector<vertex_t> medians;
    std::vector<vertex_t> randomized_medians;
    std::uniform_real_distribution<> distribution{
        0.0957, 0.9043};  // from Nagamochi's paper
    median_heuristic median_h;
    shift_heuristic& shift_h;
    median_heuristic_params params;

   public:
    /**
     * @brief construct a new probabilistic median heuristic object
     *
     * @param graph the instance
     * @param ordering vector where the ordering is stored
     * @param n_iterations number of iterations
     */
    probmedian_heuristic(instance& instance_, median_heuristic_params params =
                                                  median_heuristic_params())
        : instance_view(instance_),
          another_ordering(n_free),
          medians(n_free),
          randomized_medians(n_free),
          median_h(instance_),
          shift_h(median_h.get_shift_heuristic()),
          params(params) {
        for (vertex_t i = 0; i < n_free; ++i) another_ordering[i] = i;
        fill_medians();
    }

    // delete copy and move constructor and copy and move assignment
    probmedian_heuristic(const probmedian_heuristic& other) = delete;
    probmedian_heuristic(probmedian_heuristic&& other) = delete;
    probmedian_heuristic& operator=(probmedian_heuristic& other) = delete;
    probmedian_heuristic& operator=(const probmedian_heuristic& other) = delete;

    /**
     * @brief runs the probabilistic median heuristic solver
     * and stores the result in ordering; expects an initial solution in
     * ordering
     *
     * @param ordering a valid ordering with n_crossings number of crossings
     * @param n_crossings number of crossings of ordering
     * @return crossing_number_t number of crossings
     */
    crossing_number_t operator()(std::vector<vertex_t>& ordering,
                                 crossing_number_t n_crossings) {
        PACE_DEBUG_PRINTF("start probmedian heuristic\n");
        // try to find a better solution with probabilistic median heuristic
        for (uint8_t i = 1;
             lower_bound() < n_crossings && i <= params.n_lookahead && !tle;
             ++i) {
            const crossing_number_t candidate = generate_another_ordering();
            if (candidate < n_crossings) {
                i = 0;
                n_crossings = candidate;
                std::swap(ordering, another_ordering);
            }
            PACE_DEBUG_PRINTF("%11s=%11u, %11s=%11u\n", "i", i, "pm look",
                              params.n_lookahead);
        }
        update_ordering(ordering, n_crossings);
        PACE_DEBUG_PRINTF("end   probmedian heuristic\n");
        return n_crossings;
    }

    /**
     * @brief runs the (probabilistic) median heuristic solver and stores the
     * result in ordering
     *
     * @return crossing_number_t number of crossings
     */
    crossing_number_t operator()(std::vector<vertex_t>& ordering) {
        const crossing_number_t n_crossings = median_h(ordering);
        // try to find a better solution with probabilistic median heuristic
        return (*this)(ordering, n_crossings);
    }

    /**
     * @brief compares the randomized medians of a and b
     * we break ties by considering the medians of a and b
     * and then the parity of the degrees of a and b
     * (odd degrees are smaller)
     *
     * @return true a < b (according to medians)
     * @return false a >= b (according to medians)
     */
    inline bool compare(const vertex_t& a, const vertex_t& b) const {
        assert(n_free == randomized_medians.size());
        assert(n_free == medians.size());
        assert(a < n_free && b < n_free);
        if (randomized_medians[a] < randomized_medians[b]) {
            return true;
        } else if (randomized_medians[a] > randomized_medians[b]) {
            return false;
        } else if (medians[a] < medians[b]) {  // from here random_medians[a] ==
                                               // random_medians[b]
            return true;
        } else if (medians[a] >
                   medians[b]) {  // from here medians[a] == medians[b]
            return false;
        } else if (graph.get_degree(a) % 2 == 1 &&
                   graph.get_degree(b) % 2 == 0) {
            return true;
        } else if (graph.get_degree(a) % 2 == 0 &&
                   graph.get_degree(b) % 2 == 1) {
            return false;
        } else {
            return a < b;
        }
    }

   private:
    inline crossing_number_t generate_another_ordering() {
        fill_randomized_medians();
        sort(another_ordering.begin(), another_ordering.end(),
             [=](const vertex_t& a, const vertex_t& b) -> bool {
                 return this->compare(a, b);
             });
        return shift_h(another_ordering);
    }

    /**
     * @brief returns randomized median of `graph.get_neighbors(i)`
     */
    inline vertex_t randomized_median(const std::size_t& i) {
        const auto& neighbors = graph.get_neighbors(i);
        const std::size_t nof_neighbors = graph.get_degree(i);
        if (nof_neighbors == 0) {
            return 0;
        } else {
            const std::size_t j = static_cast<std::size_t>(
                distribution(rd_generator) * nof_neighbors);
            assert(j < nof_neighbors);
            return neighbors[j];
        }
    }

    /**
     * @brief computes randomized medians of every `graph.get_neighbors(i)`
     * using `randomized_median(i)` and stores them in `randomized_medians`
     */
    void fill_randomized_medians() {
        for (std::size_t i = 0; i < n_free; ++i) {
            randomized_medians[i] = randomized_median(i);
        }
    }

    /**
     * @brief computes medians of every `graph.get_neighbors(i)` using
     * `internal::median(...)` and stores them in `medians`
     */
    inline void fill_medians() {
        for (std::size_t i = 0; i < n_free; ++i) {
            medians[i] = median(graph.get_neighbors(i));
        }
    }
};

};  // namespace pace

#endif