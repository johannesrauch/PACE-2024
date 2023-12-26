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
    const general_bipartite_graph<T> graph;
    const T n1, n1_choose_2;
    const folded_square_matrix<T> cr_matrix;

    // lp
    glp_prob* const lp;
    std::vector<T> ordering;

   public:
    template <typename T1>
    branch_and_cut(const general_instance<T1>& instance)
        : graph(instance),
          n1(graph.get_n1()),
          n1_choose_2(n1 * (n1 - 1) / 2),
          cr_matrix(graph),
          lp(glp_create_prob()),
          ordering(n1) {
        glp_set_obj_dir(lp, GLP_MIN);

        // add n1(n1-1)/2 variables
        glp_add_cols(lp, n1_choose_2);
        int k = 1;
        for (uint32_t i = 0; i < n1; ++i) {
            for (uint32_t j = i + 1; j < n1; ++j) {
                T c_ij = cr_matrix(i, j), c_ji = cr_matrix(j, i),
                  c = c_ij - c_ji;

                if (c_ij == 0) {
                    glp_set_col_bnds(lp, k, GLP_FX, 1, 1);  // fix i < j
                } else if (c_ji == 0) {
                    glp_set_col_bnds(lp, k, GLP_FX, 0, 0);  // fix j < i
                } else {
                    // set 0 <= x_ij <= 1
                    glp_set_col_bnds(lp, k, GLP_DB, 0, 1);
                }

                glp_set_obj_coef(lp, k, (double)c);
                ++k;
            }
        }
    }

    ~branch_and_cut() { glp_delete_prob(lp); }

    inline std::size_t get_variable_index(const std::size_t i,
                                          const std::size_t j) {
        assert(i != j);
        if (i < j) {
            std::size_t offset = n1_choose_2 - (n1 - i) * (n1 - i - 1) / 2;
            return offset + j - i;
        } else {
            std::size_t offset = n1_choose_2 - (n1 - j) * (n1 - j - 1) / 2;
            return offset + i - j;
        }
    }
};

};  // namespace pace2024

#endif