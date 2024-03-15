#ifndef PACE2024_SHIFT_HEURISTIC_HPP
#define PACE2024_SHIFT_HEURISTIC_HPP

#ifndef PACE2024_CONST_SHIFT_LENGTH
#define PACE2024_CONST_SHIFT_LENGTH 128
#endif

#include <vector>

#include "bipartite_graph.hpp"
#include "crossings.hpp"
#include "debug_printf.hpp"
#include "matrix.hpp"

namespace pace2024 {

template <typename T, typename R>
class shift_heuristic {
    const folded_matrix<T> &matrix;
    std::vector<T> &ordering;
    uint32_t upper_bound;

   public:
    shift_heuristic(const folded_matrix<R> &matrix,  //
                    std::vector<T> &ordering,        //
                    const uint32_t upper_bound)
        : matrix(matrix),      //
          ordering(ordering),  //
          upper_bound(upper_bound) {
        assert(matrix.get_m() == ordering.size());
        assert(matrix.get_m() > 0);
    }

    uint32_t run(const std::size_t nof_iterations = 1024) {
        PACE2024_DEBUG_PRINTF("start shift_heuristic\n");
        bool go_on = true;
        std::size_t iteration;
        for (iteration = 0; iteration < nof_iterations && go_on; ++iteration) {
            go_on = false;
            for (std::size_t i = 0; i < ordering.size(); ++i) {
                go_on |= improve(i);
            }
        }
        PACE2024_DEBUG_PRINTF("end   shift_heuristic, iterations=%llu\n", iteration);
        return upper_bound;
    }

   private:
    /**
     * @brief shifts ordering.
     * let x[1], ..., x[n] denote the elements of ordering.
     * if i < j, then after the shift we have ..., x[i+1], ..., x[j], x[i], ...
     * if i > j, then after the shift we have ..., x[i], x[j], ..., x[i-1], ...
     */
    inline void shift(std::size_t i, std::size_t j) {
        assert(0 <= i && i < ordering.size());
        assert(0 <= j && j < ordering.size());
        assert(i != j);

        const std::ptrdiff_t delta = i < j ? 1 : -1;
        const T tmp = ordering[i];
        for (std::size_t k = i; k != j; k += delta) {
            ordering[k] = ordering[k + delta];
        }
        ordering[j] = tmp;
    }

    inline bool improve(const std::size_t &i) {
        // compute crossing difference by left shifts
        std::size_t c_old{0}, c_new{0};
        for (std::size_t j = 1; j <= PACE2024_CONST_SHIFT_LENGTH; ++j) {
            if (j > i) break;

            const std::size_t k = i - j;
            c_old += matrix(ordering[k], ordering[i]);
            c_new += matrix(ordering[i], ordering[k]);

            if (c_new < c_old) {
                shift(i, k);
                assert(c_old - c_new <= upper_bound);
                upper_bound -= c_old - c_new;
                assert(number_of_crossings(matrix, ordering) == upper_bound);
                return true;
            }
        }

        c_old = 0;
        c_new = 0;
        // compute crossing difference by right shifts
        for (std::size_t j = 1; j <= PACE2024_CONST_SHIFT_LENGTH; ++j) {
            if (i + j >= ordering.size()) break;

            const std::size_t k = i + j;
            c_old += matrix(ordering[i], ordering[k]);
            c_new += matrix(ordering[k], ordering[i]);

            if (c_new < c_old) {
                shift(i, k);
                assert(c_old - c_new <= upper_bound);
                upper_bound -= c_old - c_new;
                uint32_t ref_ub = number_of_crossings(matrix, ordering);
                assert(ref_ub == upper_bound);
                return true;
            }
        }

        return false;
    }
};

};  // namespace pace2024

#endif
