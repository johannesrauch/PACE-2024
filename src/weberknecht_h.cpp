#include <fstream>

#include "heuristic/heuristic_solver.hpp"
#include "io/output.hpp"

int main(int argc, char** argv) {
    struct sigaction action;
    memset(&action, 0, sizeof(struct sigaction));
    action.sa_handler = pace::time_limit_exceeded;
    sigaction(SIGTERM, &action, NULL);

    pace::input in;
    pace::heuristic_solver solver(in);
    std::vector<pace::vertex_t> ordering;
    solver(ordering);
    PACE_DEBUG_PRINTF("%11s=%11u\n", "tle", pace::tle);
    pace::print_output(in.get_graph().get_n_fixed(), ordering);
    return 0;
}
