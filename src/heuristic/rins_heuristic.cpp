#include "heuristic/rins_heuristic.hpp"

#include "exact/branch_and_cut.hpp"

namespace pace {

crossing_number_t rins_heuristic::operator()(  //
    highs_lp &lp, std::vector<vertex_t> &ordering) {    
    PACE_DEBUG_PRINTF("start rins\n");
    std::unique_ptr<instance> subinstance_ptr(instance_.new_rins_instance(lp));
    branch_and_cut_params params;
    params.max_nodes = 100;
    params.do_rins = false;
    branch_and_cut solver(*subinstance_ptr, params);
    const crossing_number_t n_cr = solver(ordering);
    PACE_DEBUG_PRINTF("end   rins\n");
    return shift_h(ordering, n_cr);
}

}  // namespace pace
