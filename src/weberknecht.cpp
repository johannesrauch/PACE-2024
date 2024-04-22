#include "exact/exact_solver.hpp"
#include "io/output.hpp"

int main() {
    input in;
    exact_solver solver(in);
    std::vector<vertex_t> ordering;
    solver(ordering);
    print_output(in.get_graph().get_n_fixed(), ordering);
}
