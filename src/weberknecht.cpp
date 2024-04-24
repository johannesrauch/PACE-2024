#include <fstream>

#include "exact/exact_solver.hpp"
#include "io/output.hpp"

int main(int argc, char** argv) {
    if (argc < 2) {
        pace::print_help();
        return 1;
    }

    pace::input in(argv[1]);
    pace::exact_solver solver(in);
    std::vector<pace::vertex_t> ordering;
    solver(ordering);

    if (argc < 3) {
        pace::print_output(in.get_graph().get_n_fixed(), ordering);
    } else {
        std::ofstream out_f(argv[2]);
        pace::print_output(out_f, in.get_graph().get_n_fixed(), ordering);
    }
    return 0;
}
