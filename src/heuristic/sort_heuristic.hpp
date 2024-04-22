#ifndef HEURISTIC_SORT_HEURISTIC_HPP
#define HEURISTIC_SORT_HEURISTIC_HPP

#include "heuristic/shift_heuristic.hpp"

namespace pace {

class sort_heuristic : public instance_view {

public:
sort_heuristic(instance instance_) : instance_view(instance_) {}

};

};  // namespace pace

#endif
