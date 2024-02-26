#ifndef PACE2024_BIPARTITE_GRAPH_HPP
#define PACE2024_BIPARTITE_GRAPH_HPP

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <vector>

#include "instance.hpp"

namespace pace2024 {

//
// bipartite graph
//

/**
 * @brief a bipartite graph implementation
 * used for storing the pace instances
 *
 * @tparam T vertex type
 * @tparam std::enable_if_t<std::is_unsigned<T>::value>
 */
template <typename T, class = typename std::enable_if_t<std::is_unsigned<T>::value>>
class general_bipartite_graph {
   private:
    /**
     * @brief number of vertices in the fixed partite set
     */
    std::size_t n0;

    /**
     * @brief number of vertices in the free partite set
     */
    std::size_t n1;

    /**
     * @brief number of edges
     */
    std::size_t m;

    /**
     * @brief neighbors of vertices in the free layer
     * each adjacency list is sorted
     */
    std::vector<std::vector<T>> adjacency_lists;

    /**
     * @brief used for parsing the input and storing it in general_bipartite_graph
     * friend of general_bipartite_graph
     *
     * @tparam T1 vertex type of general_bipartite_graph
     * @tparam IFSTREAM input filestream
     * @param graph
     * @param input
     */
    template <typename T1, typename IFSTREAM>
    friend void parse(general_bipartite_graph<T1> &graph, IFSTREAM &input);

   public:
    using datatype = T;

    // delete copy constructor and = operator
    general_bipartite_graph(const general_bipartite_graph &rhs) = delete;
    general_bipartite_graph &operator=(const general_bipartite_graph &rhs) = delete;
    general_bipartite_graph &operator=(general_bipartite_graph &&rhs) = delete;

    /**
     * @brief reads input from filepath and constructs general_bipartite_graph
     *
     * @param filepath
     */
    general_bipartite_graph(const std::string filepath) {
        std::ifstream input(filepath, std::ios::in);
        assert(input.good());
        parse(*this, input);
        input.close();
    }

    /**
     * @brief reads input from filestream and constructs general_bipartite_graph
     *
     * @tparam IFSTREAM
     * @param input
     */
    template <typename IFSTREAM>
    general_bipartite_graph(IFSTREAM &input) {
        parse(*this, input);
    }

    /**
     * @brief returns number of vertices in the fixed layer
     *
     * @return std::size_t
     */
    std::size_t get_n0() const { return n0; }

    /**
     * @brief returns number of vertices in the free layer
     *
     * @return std::size_t
     */
    std::size_t get_n1() const { return n1; }

    /**
     * @brief returns number of edges
     *
     * @return std::size_t
     */
    std::size_t get_m() const { return m; }

    /**
     * @brief returns a constant reference to all adjacency lists
     *
     * @return const std::vector<std::vector<T>>&
     */
    const std::vector<std::vector<T>> &get_adjacency_lists() const {
        return adjacency_lists;
    }

   private:
    /**
     * @brief sorts all adjacency lists in ascending order
     */
    void sort_adjacency_lists() {
        for (auto &adjacency_list : adjacency_lists) {
            std::sort(adjacency_list.begin(), adjacency_list.end());
        }
    }
};

using uint64_bipartite_graph = general_bipartite_graph<std::uint64_t>;
using uint32_bipartite_graph = general_bipartite_graph<std::uint32_t>;
using uint16_bipartite_graph = general_bipartite_graph<std::uint16_t>;

//
// input parse function
//

/**
 * @brief parses input from filestream and stores it in graph
 *
 * @tparam T
 * @tparam IFSTREAM
 * @param graph
 * @param input
 */
template <typename T, typename IFSTREAM>
void parse(general_bipartite_graph<T> &graph, IFSTREAM &input) {
    char type_of_line = 0;
    std::string problem_descriptor, comment;

    do {
        input >> type_of_line;
        switch (type_of_line) {
            case 'c':  // comment line
            case 'C':
                std::getline(input, comment);
                break;

            case 'p':  // parameter line
            case 'P':
                input >> problem_descriptor >> graph.n0 >> graph.n1 >> graph.m;
                break;

            default:
                std::cerr
                    << "pace2024::instance::parse(): unknown line type"
                    << std::endl;
                return;
        }
    } while (type_of_line != 'p' && type_of_line != 'P');

    graph.adjacency_lists.resize(graph.n1);
    T x, y;
    for (std::size_t i = 0; i < graph.m; ++i) {
        input >> x >> y;
        --x;
        y -= graph.n0 + 1;
        graph.adjacency_lists[y].emplace_back(x);
    }

    graph.sort_adjacency_lists();
}

};  // namespace pace2024

#endif