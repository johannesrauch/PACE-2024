#ifndef PACE2024_BRANCH_AND_CUT_HPP
#define PACE2024_BRANCH_AND_CUT_HPP

#include <glpk.h>

#include "bipartite_graph.hpp"
#include "instance.hpp"
#include "matrix.hpp"

namespace pace2024 {

template <typename T>
class branch_and_cut {
   private:
    // instance
    bipartite_graph graph;
    const uint64_t n1, n1_choose_2;
    folded_square_matrix<T> cr_matrix;

    // lp
    glp_prob* lp;

   public:
    branch_and_cut(const instance& instance)
        : graph(instance),
          n1(graph.get_n1()),
          n1_choose_2(n1 * (n1 - 1) / 2),
          cr_matrix(graph),
          lp(glp_create_prob()) {
        glp_set_obj_dir(lp, GLP_MIN);
        // add n1(n1-1)/2 variables
        glp_add_cols(lp, n1_choose_2);
        for (uint64_t j = 1; j <= n1_choose_2; ++j) {
            // set 0 <= x_ij <= 1 for all variables
            glp_set_col_bnds(lp, j, GLP_FX, 0, 1);
        }

        // set objective coefficients
        uint64_t k = 1;
        for (uint64_t i = 0; i < n1; ++i) {
            for (uint64_t j = i + 1; j < n1; ++j) {
                glp_set_obj_coef(lp, k, cr_matrix(i, j) - cr_matrix(j, i));
                ++k;
            }
        }
    }

    ~branch_and_cut() { glp_delete_prob(lp); }

    uint64_t get_variable_index(const uint64_t i, const uint64_t j) {
        assert(i != j);
        if (i < j) {
            uint64_t offset = n1_choose_2 - (n1 - i) * (n1 - i - 1) / 2;
            return offset + j - i;
        } else {
            uint64_t offset = n1_choose_2 - (n1 - j) * (n1 - j - 1) / 2;
            return offset + i - j;
        }
    }
};

};  // namespace pace2024

#endif