#include "utils/index_utils.hpp"

int main() {
    const std::size_t n = 1234;
    std::size_t ref = 0;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            assert(ref == pace::flat_index(n, i, j));
            ++ref;            
        }
    }
}
