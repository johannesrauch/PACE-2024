#include "heuristic/barycenter_heuristic.hpp"

namespace pace {

crossing_number_t barycenter_heuristic::operator()(
    std::vector<vertex_t> &ordering) {
    identity(n_free, ordering);
    std::sort(ordering.begin(), ordering.end(),
              [=](const vertex_t &a, const vertex_t &b) -> bool {
                  return this->compare(a, b);
              });
    return shift_h(ordering);
}

}  // namespace pace