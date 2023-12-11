#ifndef PACE2024_INSTANCE_HPP
#define PACE2024_INSTANCE_HPP

#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace pace2024 {

class instance {
   private:
    instance(const instance &rhs) = delete;
    instance &operator=(const instance &rhs) = delete;
    instance &operator=(instance &&rhs) = delete;

    std::string problem_descriptor;
    std::uint64_t n0,  // number of vertices in A, the fixed partite set
        n1,            // number of vertices in B
        m;             // number of edges
    std::vector<std::pair<uint64_t, uint64_t>> edges;

   public:
    instance() {}

    uint64_t get_n0() { return n0; }

    uint64_t get_n1() { return n1; }

    uint64_t get_m() { return m; }

    const std::vector<std::pair<uint64_t, uint64_t>> &get_edges() {
        return edges;
    }

    template <class ISTREAM>
    void parse(ISTREAM &input) {
        char type_of_line;
        std::string comment;
        uint64_t x, y;

        do {
            input >> type_of_line;
            switch (type_of_line) {
                case 'c':  // comment line
                case 'C':
                    std::getline(input, comment);
                    break;

                case 'p':  // parameter line
                case 'P':
                    input >> problem_descriptor >> n0 >> n1 >> m;
                    break;

                default:
                    std::cerr
                        << "pace2024::instance::parse(): unknown line type"
                        << std::endl;
                    return;
            }
        } while (type_of_line != 'p');

        edges.reserve(m);
        for (uint64_t i = 0; i < m; ++i) {
            input >> x >> y;
            edges.emplace_back(x, y);
        }
    }

    void parse() { parse(std::cin); }
};

};  // namespace pace2024

#endif
