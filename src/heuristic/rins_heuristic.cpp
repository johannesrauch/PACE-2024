#include "heuristic/rins_heuristic.hpp"

#include "exact/branch_and_cut.hpp"

namespace pace {

crossing_number_t rins_heuristic::operator()(  //
    const highs_lp &lp, std::vector<vertex_t> &ordering) {
    instance subinstance(graph);
    subinstance.update_ordering(get_ordering(), upper_bound);
    
    branch_and_cut_params params;
    params.max_nodes = 100;
    
}

}  // namespace pace
