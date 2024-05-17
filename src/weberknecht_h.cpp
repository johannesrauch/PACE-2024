#include <fstream>

#include "heuristic/heuristic_solver.hpp"
#include "io/output.hpp"

int main(int argc, char** argv) {
    pace::input in;
    pace::heuristic_solver solver(in);
    std::vector<pace::vertex_t> ordering;
    solver(ordering);
    pace::print_output(in.get_graph().get_n_fixed(), ordering);
    return 0;
}
