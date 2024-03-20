#ifndef PACE_SHIFT_HEURISTIC_HPP
#define PACE_SHIFT_HEURISTIC_HPP

#ifndef PACE_CONST_SHIFT_LENGTH
#define PACE_CONST_SHIFT_LENGTH 256
#endif

#include <list>
#include <vector>

#include "bipartite_graph.hpp"
#include "crossings.hpp"
#include "debug_printf.hpp"
#include "matrix.hpp"

namespace pace {

/**
 * @brief improves a given solution to a local minimum by shifting
 * 
 * @tparam T vertex type
 * @tparam R crossing matrix data type
 */
template <typename T, typename R>
class shift_heuristic {
    /// @brief instance
    const bipartite_graph<T> &graph;
    
    /// @brief crossing matrix
    const folded_matrix<R> &cr_matrix;

    /// @brief in-out ordering
    std::vector<T> &ordering;

    /// @brief internal ordering as list for fast shifting
    std::list<T> ordering_;

    /// @brief best upper bound (number of crossings)
    uint32_t upper_bound;

   public:
   /**
    * @brief sets references and initializes heuristic solver
    * 
    * @param graph instance
    * @param cr_matrix crossing matrix
    * @param ordering an ordering of the free layer
    * @param upper_bound number of crossings of `ordering`
    */
    shift_heuristic(const bipartite_graph<T> &graph,    //
                    const folded_matrix<R> &cr_matrix,  //
                    std::vector<T> &ordering,           //
                    const uint32_t upper_bound)
        : graph(graph),          //
          cr_matrix(cr_matrix),  //
          ordering(ordering),    //
          ordering_(ordering.begin(), ordering.end()),
          upper_bound(upper_bound) {
        assert(cr_matrix.get_m() == ordering.size());
        assert(cr_matrix.get_m() > 1);
        assert(number_of_crossings(graph, ordering) == upper_bound);
    }

    uint32_t run(const std::size_t nof_iterations = 1024) {
        PACE_DEBUG_PRINTF("start shift_heuristic\n");
        bool go_on = true;
        std::size_t iteration;
        for (iteration = 0; iteration < nof_iterations && go_on; ++iteration) {
            go_on = false;
            std::ptrdiff_t dist;
            for (typename std::list<T>::iterator it = ordering_.begin(); it != ordering_.end();) {
                const auto [go_on_, it_] = improve(it);
                go_on |= go_on_;
                it = it_;
                dist = std::distance(it, ordering_.end());
                (void)dist;
            }
        }
        std::copy(ordering_.begin(), ordering_.end(), ordering.begin());
        assert(ordering.size() == graph.get_n_free());
        assert(number_of_crossings(graph, ordering) == upper_bound);
        PACE_DEBUG_PRINTF("end   shift_heuristic, iterations=%llu\n", iteration);
        return upper_bound;
    }

   private:
    inline std::pair<bool, typename std::list<T>::iterator> improve(
        typename std::list<T>::iterator i) {
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

        // shift
        if (right_shift) ++l;
        // moves element at i in front of element at l
        // ordering_.splice(l, ordering_, i);
        ordering_.insert(l, *i);
        typename std::list<T>::iterator it = ordering_.erase(i);
        assert(ordering_.size() == graph.get_n_free());
        assert(upper_bound >= improvement);
        upper_bound -= improvement;
        assert(number_of_crossings(graph, std::vector<T>{ordering_.begin(), ordering_.end()}) ==
               upper_bound);
        return std::make_pair(true, it);
    }
};

};  // namespace pace

#endif
