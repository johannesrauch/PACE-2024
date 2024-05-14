#ifndef PACE_IO_PARSE_INPUT_HPP
#define PACE_IO_PARSE_INPUT_HPP

#include <cassert>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "model/bipartite_graph.hpp"

namespace pace {

namespace fs = std::filesystem;

/**
 * @brief creates ifstream `input` from filepath and calls `parse_input(input, graph)`
 *
 * @tparam T vertex type
 * @param input
 * @param graph in-out parameter; where instance gets stored
 */
template <typename T>
void parse_input(const fs::path filepath, general_bipartite_graph<T> &graph) {
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
template <typename CharT, typename Traits, typename T>
void parse_input(std::basic_istream<CharT, Traits> &input, general_bipartite_graph<T> &graph) {
    assert(input.good());
    graph.clear();

    std::size_t n_fixed = 0, n_free = 0, m = 0, cw = 0;
    bool param_track = false;
    char type_of_line = 0;
    std::string problem_descriptor, line;

    // parameter and comment lines
    do {
        input >> type_of_line;
        switch (type_of_line) {
            case 'c':  // comment line
            case 'C': {
                std::getline(input, line);
                break;
            }

            case 'p':  // parameter line
            case 'P': {
                std::getline(input, line);
                std::istringstream stream(line);
                stream >> problem_descriptor >> n_fixed >> n_free >> m;
                if (stream) {
                    stream >> cw;
                    param_track = true;
                }
                break;
            }

            default: {
                std::cerr << "pace::parse_input(): unknown line type" << std::endl;
                return;
            }
        }
    } while (type_of_line != 'p' && type_of_line != 'P');

    if (param_track) {
        for (std::size_t i = 0; i < n_fixed + n_free; ++i) std::getline(input, line);
    }

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

};  // namespace pace

#endif
