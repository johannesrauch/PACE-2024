#include "transitive_hull.hpp"

#include <algorithm>

#include "random.hpp"
#include "test_utils.hpp"

template <typename T>
void test_transitive_hull_with(const pace::digraph<T> &graph) {
    std::vector<std::pair<T, T>> test1, test2;

    std::clock_t start = std::clock();
    pace::transitive_hull(graph, test1);
    std::clock_t end = std::clock();
    const double t_h = pace::test::time_in_ms(start, end);

    start = std::clock();
    pace::transitive_hull_of_acyclic(graph, test2);
    end = std::clock();
    const double t_a = pace::test::time_in_ms(start, end);

    std::sort(test1.begin(), test1.end());
    std::sort(test2.begin(), test2.end());
    assert(test1 == test2);

    fmt::printf("%11u|%11.1f%11.1f\n", graph.get_n(), t_h, t_a);
    std::cout << std::flush;
}

void test_transitive_hull() {
    using T = uint16_t;
    constexpr std::size_t ns[3] = {10, 100, 1000};
    constexpr std::size_t t = 4;
    constexpr double p = 0.1;

    fmt::printf("%11s|%11s%11s\n", "n", "thull", "thull acyc");
    pace::test::print_line(34);
    for (const std::size_t &n : ns) {
        for (std::size_t i = 0; i < t; ++i) {
            std::vector<T> ordering(n);
            pace::digraph<T> graph(n);
            pace::test::shuffle(ordering);
            for (T i = 0; i < n; ++i) {
                for (T j = i + 1; j < n; ++j) {
                    if (pace::coinflip(p)) {
                        graph.add_arc(ordering[i], ordering[j]);
                    }
                }
            }

            test_transitive_hull_with(graph);
        }
    }
}

int main() {
    test_transitive_hull();
    return 0;
}