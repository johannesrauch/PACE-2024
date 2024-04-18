#ifndef PACE_EXACT_HIGHS_BASE_HPP
#define PACE_EXACT_HIGHS_BASE_HPP

#include "Highs.h"
#include "model/instance.hpp"

namespace pace {

class highs_base : public instance_view {
   protected:
    Highs solver;
    HighsStatus status{HighsStatus::kOk};

    highs_base(instance &instance_) : instance_view(instance_) {
        status = solver.setOptionValue("presolve", "off");
        assert(status == HighsStatus::kOk);
        status = solver.setOptionValue("solver", "simplex");
        assert(status == HighsStatus::kOk);
        status = solver.setOptionValue("parallel", "off");
        assert(status == HighsStatus::kOk);
        status = solver.setOptionValue("log_to_console", false);
        assert(status == HighsStatus::kOk);
    }
};

};  // namespace pace

#endif
