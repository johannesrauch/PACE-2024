#ifndef PACE2024_INPUT_HPP
#define PACE2024_INPUT_HPP

#include <cassert>
#include <cstdint>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "bipartite_graph.hpp"

namespace pace2024 {

namespace fs = std::filesystem;

/**
 * @brief creates ifstream `input` from filepath and calls `parse_input(input, graph)` 
 * 
 * @tparam T vertex type
 * @param input
 * @param graph in-out parameter; where instance gets stored
 */
template<typename T>
void parse_input(fs::path filepath, bipartite_graph<T> &graph) {
    std::ifstream input(filepath);
    parse_input(input, graph);
}

/**
 * @brief reads an instance of one-sided crossing minimization from `input` and stores it in `graph`
 * 
 * @tparam T vertex type
 * @param input 
 * @param graph in-out parameter; where instance gets stored
 */
template <typename T>
void parse_input(std::ifstream &input, bipartite_graph<T> &graph) {
    assert(input.good());
    graph.clear();

    std::size_t n_fixed = 0, n_free = 0, m = 0;
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
                input >> problem_descriptor >> n_fixed >> n_free >> m;
                break;

            default:
                std::cerr
                    << "pace2024::parse_input(): unknown line type"
                    << std::endl;
                return;
        }
    } while (type_of_line != 'p' && type_of_line != 'P');

    graph.add_fixed_vertices(n_fixed);
    graph.add_free_vertices(n_free);

    // edges
    graph.get_edges().reserve(m);
    T x, y;
    for (std::size_t i = 0; i < m; ++i) {
        input >> x >> y;
        --x;
        y -= n_fixed + 1;
        graph.add_edge(x, y);
    }
    assert(m == graph.get_m());

    // sort ascending
    graph.sort_adjacency_lists();
}

};  // namespace pace2024

#endif
