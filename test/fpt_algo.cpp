#include "fpt_algo.hpp"

#include <cstdint>

int main() {
    pace2024::general_instance<uint16_t> instance("tiny_test_set/website_20.gr");
    pace2024::fpt_algo<uint16_t> solver(instance);
    return 0;
}