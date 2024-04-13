#ifndef PACE_RANDOM_HPP
#define PACE_RANDOM_HPP

#ifndef PACE_CONST_SEED
#define PACE_CONST_SEED 42
#endif

#include <random>

namespace pace {

std::mt19937 generator(PACE_CONST_SEED);

/**
 * @brief returns true with probability p
 */
inline bool coinflip(const double p = 0.5) {
    assert(0. <= p);
    // assert(p <= 1.);
    std::uniform_real_distribution<> distribution(0., 1.);
    return distribution(generator) < p;
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

};  // namespace pace

#endif