#ifndef PACE_UTILS_RANDOMNESS_UTILS_HPP
#define PACE_UTILS_RANDOMNESS_UTILS_HPP

#ifndef PACE_CONST_SEED
#define PACE_CONST_SEED 42
#endif

#include <random>
#include "utils/vector_utils.hpp"

namespace pace {

std::mt19937 rd_generator(PACE_CONST_SEED);

/**
 * @brief returns true with probability p (if p in [0; 1])
 */
inline bool coinflip(const double p = 0.5) {
    // assert(0. <= p);
    // assert(p <= 1.);
    std::uniform_real_distribution<> distribution(0., 1.);
    return distribution(rd_generator) < p;
}

namespace test {

/**
 * @brief fills the given vector with 0, 1, ..., vec.size()-1
 * and suffles it (each permutation has some probability to appear)
 * 
 * @tparam T element type
 * @tparam A allocator
 * @param vec 
 */
template <typename T, typename A>
void shuffle(std::vector<T, A>& vec) {
    // if vec is empty, do nothing
    if (vec.size() == 0) return;

    // initialize
    for (std::size_t i = 0; i < vec.size(); ++i) {
        vec[i] = i;
    }

    // shuffle
    for (std::size_t i = 0; i < vec.size(); ++i) {
        std::uniform_int_distribution<std::size_t> distribution(i, vec.size() - 1);
        std::size_t j = distribution(rd_generator);
        T tmp = vec[i];
        vec[i] = vec[j];
        vec[j] = tmp;
    }

    assert(test::is_permutation(vec));
}

};

};  // namespace pace

#endif