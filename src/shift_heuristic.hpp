#ifndef PACE_SHIFT_HEURISTIC_HPP
#define PACE_SHIFT_HEURISTIC_HPP

#ifndef PACE_CONST_SHIFT_LENGTH
#define PACE_CONST_SHIFT_LENGTH 256
#endif

#ifndef PACE_CONST_SHIFT_ITERATIONS
#define PACE_CONST_SHIFT_ITERATIONS 4096
#endif

#include <list>
#include <vector>
#include <limits>

#include "crossings.hpp"
#include "debug_printf.hpp"
#include "instance.hpp"

namespace pace {

/**
 * @brief improves a given solution to a local minimum by shifting
 *
 * @tparam T vertex type
 * @tparam R crossing numbers type
 */
template <typename T, typename R>
class shift_heuristic {
    /// @brief instance graph
    const bipartite_graph<T> &graph;

    /// @brief crossing matrix
    const folded_matrix<R> &cr_matrix;

    /// @brief internal ordering as list for fast shifting
    std::list<T> ordering_;

    const uint32_t &lower_bound;

    /// @brief best upper bound (number of crossings)
    uint32_t upper_bound{std::numeric_limits<uint32_t>::max()};

   public:
    /**
     * @brief initializes shift heuristic for `instance`
     */
    shift_heuristic(const instance<T, R> &instance)
        : graph(instance.graph()),  //
          cr_matrix(instance.cr_matrix()),
          ordering_(graph.get_n_free()),
          lower_bound(instance.get_lower_bound()) {}

    /**
     * @brief runs the shift heuristic
     *
     * @tparam NOF_ITERATIONS as template parameter for loop unrolling
     * @param ordering in-out parameter
     */
    template <std::size_t NOF_ITERATIONS = PACE_CONST_SHIFT_ITERATIONS>
    uint32_t operator()(std::vector<T> &ordering) {
        return operator()<NOF_ITERATIONS>(ordering, number_of_crossings(graph, ordering));
    }

    /**
     * @brief runs the shift heuristic
     *
     * @tparam NOF_ITERATIONS as template parameter for loop unrolling
     * @param nof_crossings number of crossings of `ordering`
     * @param ordering in-out parameter
     */
    template <std::size_t NOF_ITERATIONS = PACE_CONST_SHIFT_ITERATIONS>
    uint32_t operator()(std::vector<T> &ordering, const uint32_t nof_crossings) {
        if (lower_bound >= nof_crossings) return nof_crossings;

        upper_bound = nof_crossings;
        assert(ordering.size() == graph.get_n_free());
        std::copy(ordering.begin(), ordering.end(), ordering_.begin());
        bool go_on = true;

        std::size_t iteration;
        for (iteration = 0; iteration < NOF_ITERATIONS && go_on; ++iteration) {
            go_on = false;
            for (typename std::list<T>::iterator it = ordering_.begin(); it != ordering_.end();) {
                const auto [go_on_, it_] = improve(it);
                go_on |= go_on_;
                it = it_;
            }
        }

        std::copy(ordering_.begin(), ordering_.end(), ordering.begin());
        assert(ordering.size() == graph.get_n_free());
        assert(number_of_crossings(graph, ordering) == upper_bound);
        return upper_bound;
    }

   private:
    /**
     * @brief tries to improve the current solution by shifts of element at i
     *
     * @return std::pair<bool, typename std::list<T>::iterator> first: true iff improvement found, second: iterator the
     * next element to consider
     */
    inline std::pair<bool, typename std::list<T>::iterator> improve(typename std::list<T>::iterator i) {
        // try left shifts
        R c_old{0}, c_new{0}, improvement = 0;
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
        assert(upper_bound >= improvement);
        upper_bound -= improvement;
        assert(number_of_crossings(graph, std::vector<T>{ordering_.begin(), ordering_.end()}) == upper_bound);
        return std::make_pair(true, it);
    }
};

};  // namespace pace

#endif
