#include <filesystem>
#include <set>

#include "instance.hpp"
#include "test_utils.hpp"

namespace fs = std::filesystem;

void test_input_with(const fs::path filepath) {
    pace::input input(filepath);
    fmt::printf("%11s%11u%11d\n", filepath.filename(), input.get_n_instances(), input.exists_trivial_instance());
    std::cout << std::flush;
}

void test_input(const fs::path dirpath) {
    fmt::printf("%s\n\n", dirpath);
    fmt::printf("%11s%11s%11s\n", "instance", "n split", "triv exist");
    pace::test::print_line(34);
    std::set<fs::path> testcases;
    for (const auto& file : fs::directory_iterator(dirpath)) {
        if (file.is_regular_file()) {
            testcases.insert(file.path());
        }
    }
    for (const auto& filepath : testcases) {
        test_input_with(filepath);
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
            test_input_with(path);
        } else {
            test_input(std::string{argv[1]} + "/instances");
        }
    }

    std::cout << "TEST::PACE::INPUT:\t\t\t\tOK" << std::endl;
    return 0;
}
