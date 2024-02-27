#ifndef PACE2024_INSTANCE_HPP
#define PACE2024_INSTANCE_HPP

#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

namespace pace2024 {

template <typename T, class = typename std::enable_if_t<std::is_unsigned<T>::value>>
class general_input {
   private:
    std::string problem_descriptor;

    std::size_t n0,  // number of vertices in A, the fixed partite set
        n1,          // number of vertices in B
        m;           // number of edges

    std::vector<std::pair<T, T>> edges;

   public:
    general_input(const general_input &rhs) = delete;
    general_input &operator=(const general_input &rhs) = delete;
    general_input &operator=(general_input &&rhs) = delete;

    general_input() { parse(); }

    general_input(std::string filename) {
        std::ifstream input(filename, std::ios::in);
        assert(input.good());
        parse(input);
        input.close();
    }

    std::size_t get_n0() const { return n0; }

    std::size_t get_n1() const { return n1; }

    std::size_t get_m() const { return m; }

    const std::vector<std::pair<T, T>> &get_edges() const {
        return edges;
    }

    template <class ISTREAM>
    void parse(ISTREAM &input) {
        char type_of_line = 0;
        std::string comment;

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
        } while (type_of_line != 'p' && type_of_line != 'P');

        edges.reserve(m);
        T x, y;
        for (uint64_t i = 0; i < m; ++i) {
            input >> x >> y;
            edges.emplace_back(x, y);
        }
    }

    void parse() { parse(std::cin); }
};

using uint64_instance = general_input<std::uint64_t>;
using uint32_instance = general_input<std::uint32_t>;
using uint16_instance = general_input<std::uint16_t>;

};  // namespace pace2024

#endif
