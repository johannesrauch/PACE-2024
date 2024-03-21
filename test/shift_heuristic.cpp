#include "shift_heuristic.hpp"

#include "barycenter_heuristic.hpp"
#include "crossings.hpp"
#include "instance.hpp"
#include "median_heuristic.hpp"
#include "test_utils.hpp"

namespace fs = std::filesystem;

template <typename T, typename R>
void test_shift_heuristic(const pace::instance<T, R> &instance, bool print_optimal) {
    std::vector<T> ordering;

    pace::barycenter_heuristic barycenter_h{instance};
    const uint32_t c_b = barycenter_h.template operator()<false>(ordering);
    const uint32_t c_sb = barycenter_h.shift_h(ordering, c_b);

    pace::median_heuristic median_h{instance};
    const uint32_t c_m = median_h.template operator()<false>(ordering);
    const uint32_t c_sm = median_h.shift_h(ordering, c_m);

    pace::probmedian_heuristic probmedian_h{instance};
    const uint32_t c_p = probmedian_h(ordering);

    int32_t optimal = -1;
    if (print_optimal) {
        optimal = pace::test::get_ref_nof_crossings(instance.filepath);
    }
    fmt::printf("%11s%11d|%11u%11u%11.4f%11u|%11u%11u%11.4f%11u|%11u%11.4f\n",
                instance.filepath.filename(), optimal,  //
                c_b, c_sb, print_optimal ? (static_cast<double>(c_sb) / optimal - 1) * 100 : -1, c_b - c_sb,
                c_m, c_sm, print_optimal ? (static_cast<double>(c_sm) / optimal - 1) * 100 : -1, c_m - c_sm,
                c_p, print_optimal ? (static_cast<double>(c_p) / optimal - 1) * 100 : -1);
    std::cout << std::flush;
}

void test_shift_heuristic_with(const fs::path dirpath, bool print_optimal) {
    fmt::printf("start test_shift_heuristic_with(%s)\n\n", dirpath);
    fmt::printf("%11s%11s|%11s%11s%11s%11s|%11s%11s%11s%11s|%11s%11s\n", "instance",
                "optimal", "barycenter", "shift", "off by %%", "delta", "median", "shift",
                "off by %%", "delta", "probmedian", "off by %%");
    pace::test::print_line(136);
    for (const auto &file : fs::directory_iterator(dirpath)) {
        if (!file.is_regular_file()) continue;

        pace::instance instance(file.path());
        test_shift_heuristic(instance, print_optimal);
    }
    fmt::printf("\nend   test_shift_heuristic_with(%s)\n\n", dirpath);
}

int main(int argc, char **argv) {
    if (argc == 1) {
        test_shift_heuristic_with("medium_test_set/instances", true);
    } else {
        test_shift_heuristic_with(std::string{argv[1]} + "/instances", false);
    }
    std::cout << "TEST::PACE::SHIFT_HEURISTIC:\t\tOK" << std::endl;
    return 0;
}
