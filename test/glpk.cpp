#include <glpk.h>

#include <iostream>

int main() {
    glp_prob *lp = glp_create_prob();
    glp_delete_prob(lp);
    std::cout << "TEST::PACE2024::GLPK:\t\t\t\tOK" << std::endl;
    return 0;
}