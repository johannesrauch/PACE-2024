#ifndef PACE2024_BIPARTITE_GRAPH_HPP
#define PACE2024_BIPARTITE_GRAPH_HPP

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace pace2024 {

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
     * @brief edges
     * first vertex is in fixed layer
     * second vertex is in free layer
     */
    std::vector<std::pair<T, T>> edges;

   public:
    using vertextype = T;

    // delete copy constructor, move constructor, copy assignment and move assignment
    general_bipartite_graph(const general_bipartite_graph &rhs) = delete;
    general_bipartite_graph(general_bipartite_graph &&rhs) = delete;
    general_bipartite_graph &operator=(const general_bipartite_graph &rhs) = delete;
    general_bipartite_graph &operator=(general_bipartite_graph &&rhs) = delete;

    /**
     * @brief creates input filestream from filepath,
     * reads from it to construct general_bipartite_graph
     *
     * @param filepath
     */
    general_bipartite_graph(const std::string filepath) {
        std::ifstream input(filepath);
        assert(input.good());
        parse_input(input);
        input.close();

        assert(n1 == adjacency_lists.size());
        assert(m == edges.size());
    }

    /**
     * @brief reads input from filestream and constructs general_bipartite_graph
     *
     * @param input
     */
    general_bipartite_graph(std::ifstream &input) {
        parse_input(input);

        assert(n1 == adjacency_lists.size());
        assert(m == edges.size());
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

    /**
     * @brief returns a reference to all edges
     *
     * @return const std::vector<std::pair<T, T>>&
     */
    std::vector<std::pair<T, T>> &get_edges() {
        return edges;
    }

    /**
     * @brief returns a constant reference to all edges
     *
     * @return const std::vector<std::pair<T, T>>&
     */
    const std::vector<std::pair<T, T>> &get_edges() const {
        return edges;
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

    /**
     * @brief used for parsing the input and storing it in general_bipartite_graph
     *
     * @param input
     */
    void parse_input(std::ifstream &input) {
        char type_of_line = 0;
        std::string problem_descriptor, comment;

        // parameter and comment lines
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
                        << "pace2024::parse_input(): unknown line type"
                        << std::endl;
                    return;
            }
        } while (type_of_line != 'p' && type_of_line != 'P');

        // edges
        adjacency_lists.resize(n1);
        edges.reserve(m);
        T x, y;
        for (std::size_t i = 0; i < m; ++i) {
            input >> x >> y;
            --x;
            y -= n0 + 1;

            adjacency_lists[y].emplace_back(x);
            edges.emplace_back(x, y);
        }

        // sort ascending
        sort_adjacency_lists();
    }
};

using uint64_bipartite_graph = general_bipartite_graph<std::uint64_t>;
using uint32_bipartite_graph = general_bipartite_graph<std::uint32_t>;
using uint16_bipartite_graph = general_bipartite_graph<std::uint16_t>;

};  // namespace pace2024

#endif