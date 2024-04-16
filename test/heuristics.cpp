#include <ctime>
#include <filesystem>
#include <fstream>

#include "barycenter_heuristic.hpp"
#include "crossings.hpp"
#include "instance.hpp"
#include "median_heuristic.hpp"
#include "oracle.hpp"
#include "printf.hpp"
#include "test_utils.hpp"
#include "vector_utils.hpp"

namespace fs = std::filesystem;

template <typename T, typename R>
uint32_t nof_dispositions(const pace::instance<T, R>& instance, std::vector<T>& ordering, const uint32_t upper_bound) {
    std::vector<T> positions;
    pace::inverse(ordering, positions);
    pace::oracle oracle(instance, upper_bound);

    const std::size_t n = instance.graph().get_n_free();
    const pace::folded_matrix<R>& cr_matrix = instance.cr_matrix();

    uint32_t nof_disp = 0;
    for (T u = 0; u < n; ++u) {
        for (T v = u + 1; v < n; ++v) {
            const double c_uv = cr_matrix(u, v);
            const double c_vu = cr_matrix(v, u);
            if (c_uv == c_vu) continue;
            pace::pattern p = oracle.based_on_crossing_numbers(c_uv, c_vu);
            if (p == pace::indeterminate) p = oracle.based_on_pattern(u, v, c_uv, c_vu);

            if (p == pace::u_before_v && positions[u] > positions[v]) ++nof_disp; 
            if (p == pace::v_before_u && positions[u] < positions[v]) ++nof_disp;
        }
    }
    return nof_disp;
}

template <typename T, typename R>
void test_heuristics(const pace::instance<T, R>& instance) {
    std::vector<T> ordering;

    // barycenter
    std::clock_t start = std::clock();
    const uint32_t c_b = pace::barycenter_heuristic{instance}(ordering);
    std::clock_t end = std::clock();
    const double t_b = pace::test::time_in_ms(start, end);
    const uint32_t d_b = nof_dispositions(instance, ordering, c_b);

    // median
    start = std::clock();
    const uint32_t c_m = pace::median_heuristic{instance}(ordering);
    end = std::clock();
    const double t_m = pace::test::time_in_ms(start, end);
    const uint32_t d_m = nof_dispositions(instance, ordering, c_m);

    // probabilistic median
    start = std::clock();
    const uint32_t c_p = pace::probmedian_heuristic{instance}(ordering);
    end = std::clock();
    const double t_p = pace::test::time_in_ms(start, end);
    assert(c_p <= c_m);
    const uint32_t d_p = nof_dispositions(instance, ordering, c_p);

    const std::string best = c_b < c_p ? "b" : (c_b == c_p ? "=" : (c_m == c_p ? "m" : "p"));
    const std::string fastest = t_b < t_m ? "b" : "m";
    const uint32_t optimal = pace::test::get_ref_n_crossings(instance.filepath);
    fmt::printf("%11s%11u%11u|%11u%11.3f%11.4f|%11u%11.3f%11.4f|%11u%11.3f%11.4f%11u|%11s%11.4f%11s|%11u%11u%11u\n",
                instance.filepath.filename(), instance.get_lower_bound(), optimal, c_b, t_b,
                (static_cast<double>(c_b) / optimal - 1) * 100, c_m, t_m,
                (static_cast<double>(c_m) / optimal - 1) * 100, c_p, t_p,
                (static_cast<double>(c_p) / optimal - 1) * 100, c_m - c_p, best,
                (static_cast<double>(std::min(c_p, c_b)) / optimal - 1) * 100, fastest,  //
                d_b, d_m, d_p);
    std::cout << std::flush;
}

void test_heuristics_with(const fs::path dirpath) {
    fmt::printf("%11s%11s%11s|%11s%11s%11s|%11s%11s%11s|%11s%11s%11s%11s|%11s%11s%11s|%11s%11s%11s\n", "instance", "lb",
                "optimal", "barycenter", "t in ms", "off by %%", "median", "t in ms", "off by %%", "probmedian",
                "t in ms", "off by %%", "delta", "best", "off by %%", "fastest", "disp b", "disp m", "disp p");
    pace::test::print_line(183);
    for (const auto& file : fs::directory_iterator(dirpath)) {
        if (!file.is_regular_file()) continue;

        pace::instance instance(file.path());
        test_heuristics(instance);
    }
    fmt::printf("\n");
}

int main(int argc, char** argv) {
    if (argc <= 1) {
        test_heuristics_with("medium_test_set/instances");
    } else {
        test_heuristics_with(std::string{argv[1]} + "/instances");
    }
    std::cout << "TEST::PACE::HEURISTICS:\t\tOK" << std::endl;
    return 0;
}