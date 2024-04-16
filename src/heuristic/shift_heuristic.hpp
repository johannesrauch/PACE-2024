#ifndef PACE_HEURISTIC_SHIFT_HEURISTIC_HPP
#define PACE_HEURISTIC_SHIFT_HEURISTIC_HPP

#ifndef PACE_CONST_SHIFT_LENGTH
#define PACE_CONST_SHIFT_LENGTH 256
#endif

#ifndef PACE_CONST_N_SHIFT_ITERATIONS
#define PACE_CONST_N_SHIFT_ITERATIONS 4096
#endif

#include <limits>
#include <list>
#include <vector>

#include "utils/crossings.hpp"
#include "debug_printf.hpp"
#include "instance.hpp"

namespace pace {

struct shift_heuristic_params {
    const std::size_t n_iterations{PACE_CONST_N_SHIFT_ITERATIONS};
    const std::size_t shift_length{PACE_CONST_SHIFT_LENGTH};
};

/**
 * @brief improves a given solution to a local minimum by shifting
 *
 * @tparam T vertex type
 * @tparam R crossing numbers type
 */
template <typename T, typename R>
class shift_heuristic : public instance_view<T, R> {
    /**
     * @brief internal ordering as list for fast shifting
     */
    std::list<T> ordering_;

    /**
     * @brief parameters
     */
    shift_heuristic_params params;

   public:
    /**
     * @brief initializes shift heuristic for `instance`
     */
    shift_heuristic(instance<T, R> &instance_,
                    shift_heuristic_params params = shift_heuristic_params())
        : instance_view<T, R>(instance_),  //
          ordering_(instance_view<T, R>::graph.get_n_free()),
          params(params) {}

    /**
     * @brief runs the shift heuristic
     *
     * @param ordering in-out parameter
     */
    uint32_t operator()(std::vector<T> &ordering) {
        return (*this)(ordering, number_of_crossings(graph, ordering));
    }

    /**
     * @brief runs the shift heuristic
     *
     * @param n_crossings number of crossings of `ordering`
     * @param ordering in-out parameter
     */
    uint32_t operator()(std::vector<T> &ordering, const uint32_t n_crossings) {
        if (lower_bound >= n_crossings) return n_crossings;

        assert(ordering.size() == graph.get_n_free());
        std::copy(ordering.begin(), ordering.end(), ordering_.begin());
        bool go_on = true;

        std::size_t iteration;
        for (iteration = 0; iteration < params.n_iterations && go_on; ++iteration) {
            go_on = false;
            for (typename std::list<T>::iterator it = ordering_.begin(); it != ordering_.end();) {
                const auto [go_on_, it_] = improve(it, n_crossings);
                go_on |= go_on_;
                it = it_;
            }
        }

        std::copy(ordering_.begin(), ordering_.end(), ordering.begin());
        assert(ordering.size() == graph.get_n_free());
        assert(number_of_crossings(graph, ordering) == n_crossings);
        upper_bound = std::min(upper_bound, n_crossings);
        
        return n_crossings;
    }

   private:
    /**
     * @brief tries to improve the current solution by shifts of element at i
     *
     * @return std::pair<bool, typename std::list<T>::iterator> first: true iff improvement found,
     * second: iterator the next element to consider
     */
    inline std::pair<bool, typename std::list<T>::iterator>  //
    improve(typename std::list<T>::iterator i, uint32_t &n_crossings) {
        // try left shifts
        uint32_t c_old{0}, c_new{0}, improvement = 0;
        typename std::list<T>::iterator k = i, l = i;
        for (std::size_t j = 1; j <= PACE_CONST_SHIFT_LENGTH; ++j) {
            if (k == ordering_.begin()) break;
            --k;

            c_old += cr_matrix(*k, *i);
            c_new += cr_matrix(*i, *k);

            if (c_new < c_old && c_old - c_new > improvement) {
                improvement = c_old - c_new;
                l = k;
            }
        }

        // try right shifts
        c_old = 0;
        c_new = 0;
        k = i;
        bool right_shift = false;
        for (std::size_t j = 1; j <= PACE_CONST_SHIFT_LENGTH; ++j) {
            ++k;
            if (k == ordering_.end()) break;

            c_old += cr_matrix(*i, *k);
            c_new += cr_matrix(*k, *i);

            if (c_new < c_old && c_old - c_new > improvement) {
                improvement = c_old - c_new;
                l = k;
                right_shift = true;
            }
        }

        if (improvement == 0) return std::make_pair(false, ++i);

        // move element at i in front of element at l
        if (right_shift) ++l;
        ordering_.insert(l, *i);
        typename std::list<T>::iterator it = ordering_.erase(i);

        assert(ordering_.size() == graph.get_n_free());
        assert(n_crossings >= improvement);
        n_crossings -= improvement;
        assert(number_of_crossings(graph, std::vector<T>{ordering_.begin(), ordering_.end()}) ==
               n_crossings);

        return std::make_pair(true, it);
    }
};

};  // namespace pace

#endif