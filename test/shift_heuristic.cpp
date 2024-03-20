#include "shift_heuristic.hpp"

#include "barycenter_heuristic.hpp"
#include "crossings.hpp"
#include "input.hpp"
#include "median_heuristic.hpp"
#include "test_utils.hpp"

namespace fs = std::filesystem;

void test_shift_heuristic(const fs::path filepath, bool print_optimal) {
    pace::bipartite_graph graph;
    pace::parse_input(filepath, graph);
    pace::folded_matrix matrix(graph.get_n_free());
    pace::fill_crossing_matrix(graph, matrix);

    std::vector<uint16_t> ordering;
    pace::barycenter_heuristic(graph, ordering).run();
    const uint32_t c_b = pace::number_of_crossings(graph, ordering);
    const uint32_t c_sb = pace::shift_heuristic(graph, matrix, ordering, c_b).run();

    pace::median_heuristic(graph, ordering).run();
    const uint32_t c_m = pace::number_of_crossings(graph, ordering);
    const uint32_t c_sm = pace::shift_heuristic(graph, matrix, ordering, c_m).run();

    pace::probmedian_heuristic heuristic(graph, ordering);
    const uint32_t c_p = heuristic.run();
    const uint32_t c_sp = pace::shift_heuristic(graph, matrix, ordering, c_p).run();

    int32_t optimal = -1;
    if (print_optimal) {
        optimal = pace::test::get_ref_nof_crossings(filepath);
    }
    fmt::printf("%11s%11d|%11u%11u%11.4f%11u|%11u%11u%11.4f%11u|%11u%11u%11.4f%11u\n",
                filepath.filename(), optimal, c_b, c_sb,
                optimal != -1 ? (static_cast<double>(c_sb) / optimal - 1) * 100 : -1, c_b - c_sb,
                c_m, c_sm, optimal != -1 ? (static_cast<double>(c_sm) / optimal - 1) * 100 : -1,
                c_m - c_sm, c_p, c_sp,
                optimal != -1 ? (static_cast<double>(c_sp) / optimal - 1) * 100 : -1, c_p - c_sp);
    std::cout << std::flush;
}

void test_shift_heuristic_with(const fs::path dirpath, bool print_optimal) {
    fmt::printf("start test_shift_heuristic_with(%s)\n\n", dirpath);
    fmt::printf("%11s%11s|%11s%11s%11s%11s|%11s%11s%11s%11s|%11s%11s%11s%11s\n", "instance",
                "optimal", "barycenter", "shift", "off by %%", "delta", "median", "shift",
                "off by %%", "delta", "probmedian", "shift", "off by %%", "delta");
    pace::test::print_line(158);
    for (const auto& file : fs::directory_iterator(dirpath)) {
        if (!file.is_regular_file()) continue;
        test_shift_heuristic(file.path(), print_optimal);
    }
    fmt::printf("\nend   test_shift_heuristic_with(%s)\n\n", dirpath);
}

int main(int argc, char** argv) {
    if (argc == 1) {
        test_shift_heuristic_with("medium_test_set/instances", true);
    } else {
        test_shift_heuristic_with("big_test_set/instances", false);
    }
    std::cout << "TEST::PACE::SHIFT_HEURISTIC:\t\tOK" << std::endl;
    return 0;
}