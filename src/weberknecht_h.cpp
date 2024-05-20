#include <fstream>

#include "heuristic/heuristic_solver.hpp"
#include "io/output.hpp"

int main(int argc, char** argv) {
    std::signal(SIGTERM, pace::timelimit::sigterm_sent);
    pace::input in;
    pace::heuristic_solver solver(in);
    std::vector<pace::vertex_t> ordering;
    solver(ordering);
    PACE_DEBUG_PRINTF("%11s=%11u\n", "tle", pace::timelimit::was_sigterm_sent());
    pace::print_output(in.get_graph().get_n_fixed(), ordering);
    return 0;
}
