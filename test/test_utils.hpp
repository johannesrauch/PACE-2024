#ifndef PACE2024_TEST_UTILS_HPP
#define PACE2024_TEST_UTILS_HPP

#include <filesystem>

namespace pace2024 {

namespace test {

template <typename R = uint32_t>
R get_ref_nof_crossings(std::filesystem::path filepath_instance) {
    auto filepath_nof_crossings =
        filepath_instance.parent_path() / "nof_crossings" / filepath_instance.filename();
    filepath_nof_crossings.replace_extension(".txt");
    // fmt::printf("%s\n", static_cast<std::string>(filepath_nof_crossings));
    std::ifstream file_nof_crossings(filepath_nof_crossings);
    assert(file_nof_crossings.good());
    R ref_nof_crossings;
    file_nof_crossings >> ref_nof_crossings;
    file_nof_crossings.close();
    return ref_nof_crossings;
}

};  // namespace test

};  // namespace pace2024

#endif