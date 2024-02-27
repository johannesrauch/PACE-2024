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

namespace test {

/**
 * @brief fills the given vector with 0, 1, ..., vec.size()-1
 * and suffles it (each permutation has some probability to appear)
 * 
 * @tparam T 
 * @param vec 
 */
template <typename T>
void shuffle(std::vector<T>& vec) {
    // if vec is empty, do nothing
    if (vec.size() == 0) return;

    // initialize
    for (std::size_t i = 0; i < vec.size(); ++i) {
        vec[i] = i;
    }

    // shuffle
    for (std::size_t i = 0; i < vec.size(); ++i) {
        std::uniform_int_distribution<std::size_t> distribution(i, vec.size() - 1);
        std::size_t j = distribution(generator);
        T tmp = vec[i];
        vec[i] = vec[j];
        vec[j] = tmp;
    }
}

};

};  // namespace pace2024

#endif