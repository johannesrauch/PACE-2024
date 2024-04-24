#ifndef PACE_IO_PIPELINE_HPP
#define PACE_IO_PIPELINE_HPP

#include <chrono>
#include <fstream>
#include <iostream>

namespace pace {

std::string generic_log_filename() {
    const auto now = std::chrono::system_clock::now();
    const std::time_t t = std::chrono::system_clock::to_time_t(now);
    return std::string(std::ctime(&t)) + ".log";
}

/**
 * @brief raii class to pipeline std::cout to a certain logfile
 */
class pipe_cout_to_file {
    std::ofstream out_f;
    std::streambuf* buf_old;

   public:
    pipe_cout_to_file(const std::string filename) : out_f(filename), buf_old(std::cout.rdbuf(out_f.rdbuf())) {}

    pipe_cout_to_file() : pipe_cout_to_file(generic_log_filename()) {}

    ~pipe_cout_to_file() { std::cout.rdbuf(buf_old); }
};

};  // namespace pace

#endif
