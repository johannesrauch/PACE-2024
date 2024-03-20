#include "debug_printf.hpp"

#include <iostream>

int main() {
    PACE_DEBUG_PRINTF("");
    std::cout << "TEST::PACE::DEBUG_PRINTF:\t\tOK" << std::endl;
    return 0;
}