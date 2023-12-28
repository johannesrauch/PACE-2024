#ifndef PACE2024_RANDOM_HPP
#define PACE2024_RANDOM_HPP

#include <random>

namespace pace2024 {

std::mt19937 generator(std::random_device{}());

/**
 * @brief returns true or false with 50% probability each
 * 
 * @return true 
 * @return false 
 */
inline bool coinflip() {
    std::uniform_real_distribution<> distribution(0., 1.);
    return distribution(generator) < 0.5;
}

};

#endif