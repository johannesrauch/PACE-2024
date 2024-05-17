#include "fmt/printf.hpp"
#include "heuristic/barycenter_heuristic.hpp"
#include "heuristic/median_heuristic.hpp"
#include "io/input.hpp"
#include "model/instance.hpp"
#include "utils/crossings_utils.hpp"
#include "utils/test_utils.hpp"
#include "utils/vector_utils.hpp"

namespace fs = std::filesystem;

using vertex_t = pace::vertex_t;
using crossing_number_t = pace::crossing_number_t;

uint32_t get_n_dispositions(pace::instance& instance, std::vector<vertex_t>& ordering) {
    std::vector<vertex_t> positions;
    pace::inverse(ordering, positions);

    const std::size_t n_free = instance.n_free;
    const pace::crossing_matrix& cr_matrix = instance.get_cr_matrix();

    uint32_t nof_disp = 0;
    for (vertex_t u = 0u; u + 1u < n_free; ++u) {
        for (vertex_t v = u + 1; v < n_free; ++v) {
            const double c_uv = cr_matrix(u, v);
            const double c_vu = cr_matrix(v, u);
            if (c_uv == c_vu) continue;
            pace::pattern p = instance.based_on_crossing_numbers(c_uv, c_vu);
            if (p == pace::indeterminate) p = instance.based_on_pattern(u, v, c_uv, c_vu);

            if (p == pace::u_before_v && positions[u] > positions[v]) ++nof_disp;
            if (p == pace::v_before_u && positions[u] < positions[v]) ++nof_disp;
        }
    }
    return nof_disp;
}

void test_heuristics(const fs::path filepath) {
    pace::input input(filepath);
    pace::instance instance(input.get_graph());
    std::vector<vertex_t> ordering;

    // barycenter
    std::clock_t start = std::clock();
    const crossing_number_t c_b = pace::barycenter_heuristic{instance}(ordering);
    std::clock_t end = std::clock();
    const double t_b = pace::test::time_in_ms(start, end);
    const uint32_t d_b = get_n_dispositions(instance, ordering);

    // median
    start = std::clock();
    const crossing_number_t c_m = pace::median_heuristic{instance}(ordering);
    end = std::clock();
    const double t_m = pace::test::time_in_ms(start, end);
    const uint32_t d_m = get_n_dispositions(instance, ordering);

    // probabilistic median
    start = std::clock();
    const crossing_number_t c_p = pace::probmedian_heuristic{instance}(ordering, c_m);
    end = std::clock();
    const double t_p = pace::test::time_in_ms(start, end);
    assert(c_p <= c_m);
    const uint32_t d_p = get_n_dispositions(instance, ordering);

    const std::string best = c_b < c_p ? "b" : (c_b == c_p ? "=" : (c_m == c_p ? "m" : "p"));
    const std::string fastest = t_b < t_m ? "b" : "m";
    const crossing_number_t optimal = pace::test::get_ref_n_crossings(filepath);
    fmt::printf("%11s%11u%11u|%11u%11.3f%11.4f|%11u%11.3f%11.4f|%11u%11.3f%11.4f%11u|%11s%11.4f%11s|%11u%11u%11u\n",
                filepath.filename(), instance.get_lower_bound(), optimal,  //
                c_b, t_b, (static_cast<double>(c_b) / optimal - 1) * 100,  //
                c_m, t_m, (static_cast<double>(c_m) / optimal - 1) * 100,  //
                c_p, t_p, (static_cast<double>(c_p) / optimal - 1) * 100, c_m - c_p, best,
                (static_cast<double>(std::min(c_p, c_b)) / optimal - 1) * 100, fastest, d_b, d_m, d_p);
    std::cout << std::flush;
}

void test_heuristics_with(const fs::path dirpath) {
    fmt::printf("%11s%11s%11s|%11s%11s%11s|%11s%11s%11s|%11s%11s%11s%11s|%11s%11s%11s|%11s%11s%11s\n",  //
                "instance", "lb", "optimal",                                                            //
                "barycenter", "t in ms", "off by %%",                                                   //
                "median", "t in ms", "off by %%",                                                       //
                "probmedian", "t in ms", "off by %%", "delta",                                          //
                "best", "off by %%", "fastest",                                                         //
                "disp b", "disp m", "disp p");
    pace::test::print_line(215);
    for (const auto& file : fs::directory_iterator(dirpath)) {
        if (!file.is_regular_file()) continue;
        test_heuristics(file.path());
    }
    fmt::printf("\n");
}

int main(int argc, char** argv) {
    if (argc <= 1) {
        test_heuristics_with("medium/instances");
    } else {
        const fs::path path = argv[1];
        if (path.extension() == ".gr") {
            test_heuristics(path);
        } else {
            test_heuristics_with(std::string{argv[1]} + "/instances");
        }
    }
    return 0;
}