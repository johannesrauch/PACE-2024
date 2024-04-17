#include "io/input.hpp"

#include <set>

#include "model/instance.hpp"
#include "utils/crossings_utils.hpp"
#include "utils/test_utils.hpp"
#include "utils/vector_utils.hpp"
#include "utils/randomness_utils.hpp"

namespace fs = std::filesystem;

using vertex_t = pace::vertex_t;
using crossing_number_t = pace::crossing_number_t;

void test_input_subgraphs(pace::input &input) {
    const std::size_t n_subgraphs = input.get_n_subgraphs();
    if (n_subgraphs < 2) return;
    std::vector<vertex_t> ordering(input.get_graph().get_n_free());

    crossing_number_t n_crossings{0};
    for (std::size_t i = 0; i < n_subgraphs; ++i) {
        const pace::bipartite_graph &subgraph = input.get_subgraph(i);
        pace::instance subinstance(subgraph);

        std::vector<vertex_t> subordering(subinstance.n_free);
        pace::test::shuffle(subordering);

        if (subinstance.n_free > 1) n_crossings += pace::number_of_crossings(subgraph, subordering);
        input.lift_ordering(i, subordering, ordering);
    }

    assert(pace::test::is_permutation(ordering));
    const crossing_number_t ref_n_crossings = pace::number_of_crossings(input.get_graph(), ordering);
    (void) ref_n_crossings;
    assert(n_crossings == ref_n_crossings);
}

void test_input(const fs::path filepath) {
    pace::input input(filepath);
    input.try_split();
    test_input_subgraphs(input);

    fmt::printf("%11s%11u%11s| ", filepath.filename(), input.get_n_subgraphs(),
                input.is_first_graph_empty() ? "true" : "false");
    for (std::size_t i = 0; i < input.get_n_subgraphs(); ++i) {
        fmt::printf("%d, ", input.get_subgraph(i).get_n_free());
    }
    fmt::printf("\n");
    std::cout << std::flush;
}

void test_input_with(const fs::path dirpath) {
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
        test_input(filepath);
    }
    fmt::printf("\n");
}

int main(int argc, char** argv) {
    if (argc <= 1) {
        test_input_with("tiny_test_set/instances");
        test_input_with("my_tests/instances");
    } else {
        const fs::path path = argv[1];
        if (path.extension() == ".gr") {
            test_input(path);
        } else {
            test_input_with(std::string{argv[1]} + "/instances");
        }
    }

    std::cout << "TEST::PACE::INPUT:\t\t\t\tOK" << std::endl;
    return 0;
}
