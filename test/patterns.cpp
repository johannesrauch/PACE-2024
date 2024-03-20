#include "patterns.hpp"

#include "bipartite_graph.hpp"
#include "input.hpp"
#include "matrix.hpp"
#include "test_utils.hpp"
#include "vector_utils.hpp"

namespace fs = std::filesystem;

using T = uint16_t;
using R = uint32_t;

void test_pattern_analyzer(const fs::path filepath, const bool test) {
    pace::bipartite_graph graph;
    pace::parse_input(filepath, graph);
    pace::folded_matrix cr_matrix(graph);
    pace::pattern_analyzer pattern(graph, cr_matrix);
    std::vector<T> ordering;
    std::vector<T> positions;
    if (test) {
        pace::test::get_ref_ordering(  //
            filepath, graph.get_n_fixed(), graph.get_n_free(), ordering);
        pace::inverse(ordering, positions);
    }

    std::vector<std::pair<T, T>> new_fixings;
    for (T u = 0; u < graph.get_n_free(); ++u) {
        for (T v = u + 1; v < graph.get_n_free(); ++v) {
            pace::pattern p = pattern.analyze(u, v);
            if (test) {
                if (p == pace::pattern::u_before_v) {
                    assert(positions[u] < positions[v]);
                } else if (p == pace::pattern::v_before_u) {
                    assert(positions[u] > positions[v]);
                }
            }
            if (p != pace::pattern::indeterminate &&  //
                cr_matrix(u, v) > 0 &&                    //
                cr_matrix(v, u) > 0) {
                new_fixings.emplace_back(u, v);
            }
        }
    }
    fmt::printf("%11s%11u%11s", filepath.filename(), new_fixings.size(), "");
    for (std::size_t i = 0; i < 3 && i < new_fixings.size(); ++i) {
        const auto& [u, v] = new_fixings[i];
        fmt::printf("(d_u=%2u,d_v=%2u,c_uv=%3u,c_vu=%3u) ",            //
                    graph.degree_of_free(u), graph.degree_of_free(v),  //
                    cr_matrix(u, v), cr_matrix(v, u));
    }
    fmt::printf("\n");
}

void test_pattern_analyzer_with(const fs::path dirpath, const bool test) {
    fmt::printf("start test_pattern_analyzer_with(%s)\n\n", dirpath);
    fmt::printf("%11s%11s\n", "instance", "new fix");
    pace::test::print_line(23);
    for (const auto& file : fs::directory_iterator(dirpath)) {
        test_pattern_analyzer(file.path(), test);
    }
    fmt::printf("\nend test_pattern_analyzer_with(%s)\n\n", dirpath);
}

int main() {
    test_pattern_analyzer_with("medium_test_set/instances", true);
    // test_pattern_analyzer_with("big_test_set/instances", false);
    std::cout << "TEST::PACE::PATTERNS:\t\t\t\t\tOK" << std::endl;
    return 0;
}
