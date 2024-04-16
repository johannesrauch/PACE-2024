#include "io/input.hpp"

#include <filesystem>
#include <set>

#include "instance.hpp"
#include "test_utils.hpp"
#include "vector_utils.hpp"

namespace fs = std::filesystem;

template <typename T, typename R>
void test_input_with(pace::input<T, R>& input) {
    input.try_split();
    fmt::printf("%11s%11u%11s| ", input.filepath.filename(), input.get_n_subgraphs(),
                input.is_first_graph_empty() ? "true" : "false");

    for (std::size_t i = 0; i < input.get_n_subgraphs(); ++i) {
        fmt::printf("%d, ", input.get_subgraph(i).get_n_free());
    }
    fmt::printf("\n");
    std::cout << std::flush;
}

void test_input(const fs::path dirpath) {
    fmt::printf("%s\n\n", dirpath);
    fmt::printf("%11s%11s%11s|%11s\n", "instance", "n subgr", "first triv", "ns free");
    pace::test::print_line(46);
    std::set<fs::path> testcases;
    for (const auto& file : fs::directory_iterator(dirpath)) {
        if (file.is_regular_file()) {
            testcases.insert(file.path());
        }
    }
    for (const auto& filepath : testcases) {
        pace::input input(filepath);
        test_input_with(input);
    }
    fmt::printf("\n");
}

int main(int argc, char** argv) {
    if (argc <= 1) {
        test_input("tiny_test_set/instances");
        test_input("my_tests/instances");
    } else {
        const fs::path path = argv[1];
        if (path.extension() == ".gr") {
            pace::input input(path);
            test_input_with(input);
        } else {
            test_input(std::string{argv[1]} + "/instances");
        }
    }

    std::cout << "TEST::PACE::INPUT:\t\t\t\tOK" << std::endl;
    return 0;
}
