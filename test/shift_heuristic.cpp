#include "shift_heuristic.hpp"

#include "barycenter_heuristic.hpp"
#include "crossings.hpp"
#include "input.hpp"
#include "median_heuristic.hpp"
#include "test_utils.hpp"

namespace fs = std::filesystem;

void test_shift_heuristic(const fs::path filepath) {
    pace2024::bipartite_graph graph;
    pace2024::parse_input(filepath, graph);
    pace2024::folded_matrix matrix(graph.get_n_free());
    pace2024::fill_crossing_matrix(graph, matrix);

    std::vector<uint16_t> ordering;
    pace2024::barycenter_heuristic(graph, ordering).run();
    const uint32_t c_b = pace2024::number_of_crossings(matrix, ordering);
    const uint32_t c_sb = pace2024::shift_heuristic(matrix, ordering, c_b).run();

    pace2024::median_heuristic(graph, ordering).run();
    const uint32_t c_m = pace2024::number_of_crossings(matrix, ordering);
    const uint32_t c_sm = pace2024::shift_heuristic(matrix, ordering, c_m).run();

    pace2024::probabilistic_median_heuristic heuristic(graph, ordering);
    const uint32_t c_p = heuristic.run();
    const uint32_t c_sp = pace2024::shift_heuristic(matrix, ordering, c_p).run();

    const uint32_t optimal = pace2024::test::get_ref_nof_crossings(filepath);
    fmt::printf("%11s%11u|%11u%11u%11.4f%11u|%11u%11u%11.4f%11u|%11u%11u%11.4f%11u\n", 
                filepath.filename(), optimal, 
                c_b, c_sb, (static_cast<double>(c_sb) / optimal - 1) * 100, c_b - c_sb,
                c_m, c_sm, (static_cast<double>(c_sm) / optimal - 1) * 100, c_m - c_sm, 
                c_p, c_sp, (static_cast<double>(c_sp) / optimal - 1) * 100, c_p - c_sp);
    std::cout << std::flush;
}

int main() {
    fmt::printf("%11s%11s|%11s%11s%11s%11s|%11s%11s%11s%11s|%11s%11s%11s%11s\n", 
                "instance", "optimal",
                "barycenter", "shift", "off by %%", "delta",
                "median", "shift", "off by %%", "delta", 
                "probmedian", "shift", "off by %%", "delta");
    pace2024::test::print_line(158);
    for (const auto& file : fs::directory_iterator("medium_test_set/instances")) {
        if (!file.is_regular_file()) continue;
        test_shift_heuristic(file.path());
    }
    return 0;
}