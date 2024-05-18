#ifndef PACE_HEURISTIC_SHIFT_HEURISTIC_HPP
#define PACE_HEURISTIC_SHIFT_HEURISTIC_HPP

#include <list>

#include "log/debug_printf.hpp"
#include "model/instance.hpp"
#include "utils/crossings_utils.hpp"
#include "utils/tle.hpp"

namespace pace {

struct shift_heuristic_params {
    const std::size_t n_iterations{4096};
    const std::size_t shift_length{512};
};

/**
 * @brief improves a given solution to a local minimum by shifting
 */
class shift_heuristic : public instance_view {
    /**
     * @brief internal ordering as list for fast shifting
     */
    std::list<vertex_t> ordering_;

    /**
     * @brief parameters
     */
    shift_heuristic_params params;

   public:
    /**
     * @brief initializes shift heuristic for `instance`
     */
    shift_heuristic(instance &instance_,
                    shift_heuristic_params params = shift_heuristic_params())
        : instance_view(instance_),  //
          ordering_(graph.get_n_free()),
          params(params) {}

    /**
     * @brief runs the shift heuristic
     *
     * @param ordering in-out parameter
     */
    crossing_number_t operator()(std::vector<vertex_t> &ordering) {
        return (*this)(ordering, number_of_crossings(graph, ordering));
    }

    /**
     * @brief runs the shift heuristic
     *
     * @param n_crossings number of crossings of 'ordering'
     * @param ordering in-out parameter
     */
    crossing_number_t operator()(  //
        std::vector<vertex_t> &ordering, crossing_number_t n_crossings);

   private:
    /**
     * @brief tries to improve the current solution by shifts of element at i
     *
     * @return std::pair<bool, typename std::list<vertex_t>::iterator> first:
     * true iff improvement found, second: iterator the next element to consider
     */
    std::pair<bool, typename std::list<vertex_t>::iterator>  //
    improve(typename std::list<vertex_t>::iterator i, uint32_t &n_crossings);
};

};  // namespace pace

#endif
