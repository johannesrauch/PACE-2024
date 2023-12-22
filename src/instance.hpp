#ifndef PACE2024_INSTANCE_HPP
#define PACE2024_INSTANCE_HPP

#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

namespace pace2024 {

template<typename T, class = typename std::enable_if_t<std::is_unsigned<T>::value>>
class general_instance {
   private:
    std::string problem_descriptor;
    T n0,  // number of vertices in A, the fixed partite set
        n1,            // number of vertices in B
        m;             // number of edges
    std::vector<std::pair<T, T>> edges;

   public:
    general_instance(const general_instance &rhs) = delete;
    general_instance &operator=(const general_instance &rhs) = delete;
    general_instance &operator=(general_instance &&rhs) = delete;

    general_instance() { parse(); }

    general_instance(std::string filename) {
        std::ifstream input(filename, std::ios::in);
        parse(input);
        input.close();
    }

    T get_n0() const { return n0; }

    T get_n1() const { return n1; }

    T get_m() const { return m; }

    const std::vector<std::pair<T, T>> &get_edges() const {
        return edges;
    }

    template <class ISTREAM>
    void parse(ISTREAM &input) {
        char type_of_line;
        std::string comment;
        T x, y;

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
        for (uint64_t i = 0; i < m; ++i) {
            input >> x >> y;
            edges.emplace_back(x, y);
        }
    }

    void parse() { parse(std::cin); }
};

using uint64_instance = general_instance<std::uint64_t>;
using uint32_instance = general_instance<std::uint32_t>;

};  // namespace pace2024

#endif
