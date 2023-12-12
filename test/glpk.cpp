#include <glpk.h>

#include <iostream>

int main() {
    glp_prob *lp = glp_create_prob();
    glp_delete_prob(lp);
    std::cout << "TEST::PACE2024::GLPK: OKAY" << std::endl;
    return 0;
}