#ifndef PACE2024_CROSSING_NUMBER_MATRIX_HPP
#define PACE2024_CROSSING_NUMBER_MATRIX_HPP

#include <cstdint>

#include "bipartite_graph.hpp"
#include "printf.hpp"

namespace pace2024 {

class crossing_number_matrix {
   private:
    uint64_t n1;
    uint64_t *data;  // row-major

    // adapted from Dujmovic and Whitesides
    void fill(bipartite_graph &graph) {
        graph.sort_adjacency_lists();
        auto adjacency_lists = graph.get_adjacency_lists();

        for (uint64_t v = 0; v < n1; ++v) {
            if (adjacency_lists[v].size() == 0) continue;
            uint64_t rv = adjacency_lists[v][adjacency_lists[v].size() - 1];

            for (uint64_t w = 0; w < n1; ++w) {
                if (v == w) continue;

                for (uint64_t wp : adjacency_lists[w]) {
                    if (wp >= rv) break;

                    // returns it to the first element such that *it >= wp + 1
                    // using binary search
                    auto it =
                        std::lower_bound(adjacency_lists[v].begin(),
                                         adjacency_lists[v].end(), wp + 1);
                    if (it != adjacency_lists[v].end()) assert(*it >= wp + 1);
                    if (it != adjacency_lists[v].begin())
                        assert(*(it - 1) < wp + 1);
                    (*this)(v, w) +=
                        std::distance(it, adjacency_lists[v].end());
                }
            }
        }
    }

   public:
    crossing_number_matrix(const crossing_number_matrix &rhs) = delete;
    crossing_number_matrix &operator=(const crossing_number_matrix &rhs) =
        delete;
    crossing_number_matrix &operator=(crossing_number_matrix &&rhs) = delete;

    crossing_number_matrix(bipartite_graph &graph)
        : n1(graph.get_n1()), data(new uint64_t[n1 * n1]()) {
        // the additional () above initializes memory to 0
        fill(graph);
    }

    crossing_number_matrix(const uint64_t n1)
        : n1(n1),
          data(new uint64_t[n1 * n1]()) {
    }  // the additional () initializes memory to 0

    ~crossing_number_matrix() { delete[] data; }

    uint64_t get_n1() const { return n1; }

    uint64_t &operator()(uint64_t i, uint64_t j) { return data[i * n1 + j]; }

    const uint64_t &operator()(uint64_t i, uint64_t j) const {
        return data[i * n1 + j];
    }
};

void fill_naivly(bipartite_graph &graph, crossing_number_matrix &matrix) {
    assert(graph.get_n1() == matrix.get_n1());
    uint64_t n1 = matrix.get_n1();
    auto adjacency_lists = graph.get_adjacency_lists();

    for (uint64_t i = 0; i < n1 - 1; ++i) {
        for (uint64_t j = i + 1; j < n1; ++j) {
            for (uint64_t x : adjacency_lists[i]) {
                for (uint64_t y : adjacency_lists[j]) {
                    if (x < y)
                        ++matrix(j, i);
                    else if (x > y)
                        ++matrix(i, j);
                }
            }
        }
    }
}

void print_crossing_number_matrix(const crossing_number_matrix &matrix) {
    uint64_t n1 = matrix.get_n1();
    for (uint64_t i = 0; i < n1; ++i) {
        for (uint64_t j = 0; j < n1; ++j) {
            // Bossert's printf does not care for the right wildcard <3
            fmt::printf("%5s", matrix(i, j));
        }
        fmt::printf("\n");
    }
}

}  // namespace pace2024

#endif