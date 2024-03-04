#include "coin-or/ClpSimplex.hpp"

int main() {
    ClpSimplex model;
    model.resize(42, 42);
    std::cout << "TEST::PACE2024::CLP:\t\t\t\tOK" << std::endl;
    return 0;
}