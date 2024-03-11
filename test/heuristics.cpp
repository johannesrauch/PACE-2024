#include <ctime>
#include <filesystem>
#include <fstream>

#include "barycenter_heuristic.hpp"
#include "bipartite_graph.hpp"
#include "crossing_number.hpp"
#include "input.hpp"
#include "median_heuristic.hpp"
#include "printf.hpp"
#include "test_utils.hpp"

namespace fs = std::filesystem;

using T = uint16_t;
using R = uint32_t;

double time_in_ms(const std::clock_t start, const std::clock_t end) {
    return 1000.0 * (end - start) / CLOCKS_PER_SEC;
}

void compare_heuristics_on_instance(fs::path filepath) {
    // input
    std::ifstream input(filepath);
    assert(input.good());
    pace2024::bipartite_graph<T> graph;
    pace2024::parse_input(filepath, graph);
    std::vector<T> ordering;

    // barycenter
    std::clock_t start = std::clock();
    pace2024::barycenter_heuristic<T>(graph, ordering).run();
    std::clock_t end = std::clock();
    const R c_b = pace2024::number_of_crossings<T, R>(graph, ordering);
    const double t_b = time_in_ms(start, end);

    // median
    start = std::clock();
    pace2024::probabilistic_median_heuristic<T, R> heuristic(graph, ordering);
    end = std::clock();
    const R c_m = heuristic.get_best();
    const double t_m = time_in_ms(start, end);

    // probabilistic median
    start = std::clock();
    const R c_p = heuristic.run(1000);
    end = std::clock();
    const double t_p = time_in_ms(start, end);
    assert(c_p <= c_m);

    const std::string best = c_b < c_p ? "b" : (c_b == c_p ? "=" : (c_m == c_p ? "m" : "p"));
    const std::string fastest = t_b < t_m ? "b" : "m";
    const R optimal = pace2024::test::get_ref_nof_crossings<R>(filepath);
    fmt::printf("%11u%11s|%11u%11.3f%11.1f|%11u%11.3f%11.1f|%11u%11.3f%11.1f%11u|%11s%11s\n",
                filepath.filename(), optimal,
                c_b, t_b, (static_cast<double>(c_b) / optimal - 1) * 100,
                c_m, t_m, (static_cast<double>(c_m) / optimal - 1) * 100,
                c_p, t_p, (static_cast<double>(c_p) / optimal - 1) * 100, c_m - c_p,
                best, fastest);
}

void compare_heuristics_on_medium_test_set() {
    fmt::printf("%11s%11s|%11s%11s%11s|%11s%11s%11s|%11s%11s%11s%11s|%11s%11s\n",
                "instance", "optimal",
                "barycenter", "t in ms", "off by %%",
                "median", "t in ms", "off by %%",
                "probmedian", "t in ms", "off by %%", "delta",
                "best", "fastest");
    fmt::printf("---------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
    for (const auto& file : fs::directory_iterator("medium_test_set")) {
        if (!file.is_regular_file()) continue;
        compare_heuristics_on_instance(file.path());
    }
}

int main() {
    compare_heuristics_on_medium_test_set();
    return 0;
}