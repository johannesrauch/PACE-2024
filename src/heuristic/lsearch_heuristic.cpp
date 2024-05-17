#include "heuristic/lsearch_heuristic.hpp"

#include "exact/branch_and_cut.hpp"

namespace pace {

crossing_number_t lsearch_heuristic::operator()(
    std::vector<vertex_t> &ordering) {
    PACE_DEBUG_PRINTF("start lsearch\n");
    std::unique_ptr<instance> subinstance_ptr(
        instance_.new_lsearch_instance(lsearch_width));
    branch_and_cut_params params;
    params.max_nodes = 50;
    params.do_rins = false;
    params.do_uninformed_h = false;
    params.do_uninformed_rins = false;
    branch_and_cut solver(*subinstance_ptr, params);
    crossing_number_t n_cr = solver(ordering);
    n_cr = shift_h(ordering, n_cr);
    PACE_DEBUG_PRINTF("end   lsearch, upper_bound=%u\n", upper_bound);
    return n_cr;
}

}  // namespace pace
