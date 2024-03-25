#include "oracle.hpp"

#include <cassert>

#include "crossings.hpp"
#include "heuristics.hpp"
#include "instance.hpp"
#include "test_utils.hpp"
#include "vector_utils.hpp"

namespace fs = std::filesystem;

template <typename T>
void insert_based_on_pattern(
    std::vector<std::pair<T, T>>& vec,
    const T& u,
    const T& v,
    pace::pattern& p) {
    assert(p != pace::indeterminate);
    if (p == pace::u_before_v) {
        vec.emplace_back(u, v);
    } else {
        vec.emplace_back(v, u);
    }
}

template <typename T, typename R>
void test_position_oracle(const pace::instance<T, R>& instance, const bool do_test) {
    const pace::bipartite_graph<T>& graph(instance.graph());
    const pace::folded_matrix<R>& cr_matrix(instance.cr_matrix());
    std::vector<T> ordering;
    std::vector<T> positions;
    const uint32_t lb = pace::get_lower_bound(cr_matrix);
    const uint32_t ub = pace::heuristics(instance, ordering);
    assert(lb <= ub);
    pace::position_oracle oracle(instance, lb, ub);

    std::vector<std::pair<T, T>> settled_pairs;
    uint32_t c = 0, b = 0, p = 0;
    const std::size_t n = graph.get_n_free();
    for (T u = 0; u < n; ++u) {
        for (T v = u + 1; v < n; ++v) {
            const R& c_uv = cr_matrix(u, v);
            const R& c_vu = cr_matrix(v, u);
            if (c_uv == c_vu) continue;

            // c_uv != c_vu from her
            pace::pattern pattern = oracle.based_on_crossing_numbers(c_uv, c_vu);
            if (pattern != pace::indeterminate) {
                ++c;
                insert_based_on_pattern(settled_pairs, u, v, pattern);
                continue;
            }

            pattern = oracle.based_on_bounds(c_uv, c_vu);
            if (pattern != pace::indeterminate) {
                ++b;
                insert_based_on_pattern(settled_pairs, u, v, pattern);
                continue;
            }

            pattern = oracle.based_on_pattern(u, v, c_uv, c_vu);
            if (pattern != pace::indeterminate) {
                ++p;
                insert_based_on_pattern(settled_pairs, u, v, pattern);
                continue;
            }
        }
    }

    if (do_test) {
        pace::test::get_ref_ordering(  //
            instance.filepath, graph.get_n_fixed(), graph.get_n_free(), ordering);
        pace::inverse(ordering, positions);
        for (const auto& [u, v] : settled_pairs) {
            assert(positions[u] < positions[v]);
        }
    }

    std::size_t n_choose_2 = n * (n - 1) / 2;
    fmt::printf("%11s|%11u%11u%11.1f|%11u%11u%11u|%11u%11u\n",
                instance.filepath.filename(),
                n_choose_2, settled_pairs.size(), static_cast<double>(settled_pairs.size()) / n_choose_2 * 100,
                c, b, p, lb, ub);
    std::cout << std::flush;
}

void test_pattern_analyzer_with(const fs::path dirpath, const bool do_test) {
    fmt::printf("start test_pattern_analyzer_with(%s)\n\n", dirpath);
    fmt::printf("%11s|%11s%11s%11s|%11s%11s%11s|%11s%11s\n",
                "instance",
                "vars", "settled", "%%",
                "crossing n", "bounds", "pattern",
                "lb", "ub");
    pace::test::print_line(103);
    for (const auto& file : fs::directory_iterator(dirpath)) {
        pace::instance instance(file.path());
        test_position_oracle(instance, do_test);
    }
    fmt::printf("\nend test_pattern_analyzer_with(%s)\n\n", dirpath);
}

/**
 * @brief tests oracle.hpp
 */
int main(int argc, char** argv) {
    if (argc <= 1) {
        test_pattern_analyzer_with("medium_test_set/instances", true);
    } else {
        test_pattern_analyzer_with(std::string{argv[1]} + "/instances", false);
    }
    std::cout << "TEST::PACE::PATTERNS:\t\t\t\t\tOK" << std::endl;
    return 0;
}
