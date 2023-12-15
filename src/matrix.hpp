#ifndef PACE2024_MATRIX_HPP
#define PACE2024_MATRIX_HPP

#include <cstdint>

#include "bipartite_graph.hpp"
#include "printf.hpp"

namespace pace2024 {

class matrix;

void fill(bipartite_graph &graph, matrix &matrix, bool clear);
void fill2(bipartite_graph &graph, matrix &cr_matrix, bool clear);

class matrix {
   private:
    uint64_t m, n;   // number of rows and columns, respectively
    uint64_t *data;  // row-major

   public:
    matrix(const matrix &rhs) = delete;
    matrix &operator=(const matrix &rhs) = delete;
    matrix &operator=(matrix &&rhs) = delete;

    matrix(bipartite_graph &graph)
        : m(graph.get_n1()),  //
          n(graph.get_n1()),
          data(new uint64_t[m * n]()) {
        // the additional () above initializes memory to 0
        assert(graph.get_n1() > 0);
        fill2(graph, *this, false);
    }

    matrix(const uint64_t m, const uint64_t n)
        : m(m), n(n), data(new uint64_t[m * n]()) {
        assert(m > 0);
        assert(n > 0);
    }  // the additional () initializes memory to 0

    ~matrix() { delete[] data; }

    uint64_t get_m() const { return m; }

    uint64_t get_n() const { return n; }

    uint64_t &operator()(uint64_t i, uint64_t j) { return data[i * n + j]; }

    const uint64_t &operator()(uint64_t i, uint64_t j) const {
        return data[i * n + j];
    }
};

// adapted from Dujmovic and Whitesides
void fill(bipartite_graph &graph, matrix &cr_matrix, bool clear = true) {
    assert(graph.get_n1() == cr_matrix.get_m());
    graph.sort_adjacency_lists();
    uint64_t n1 = graph.get_n1();
    auto adjacency_lists = graph.get_adjacency_lists();
    if (clear) memset(&cr_matrix(0, 0), 0, sizeof(uint64_t) * n1 * n1);

    for (uint64_t v = 0; v < n1; ++v) {
        if (adjacency_lists[v].size() == 0) continue;
        uint64_t rv = adjacency_lists[v][adjacency_lists[v].size() - 1];

        for (uint64_t w = 0; w < n1; ++w) {
            if (v == w) continue;

            for (uint64_t wp : adjacency_lists[w]) {
                if (wp >= rv) break;

                // returns it to the first element such that *it >= wp + 1
                // using binary search
                auto it = std::lower_bound(adjacency_lists[v].begin(),
                                           adjacency_lists[v].end(), wp + 1);
                if (it != adjacency_lists[v].end()) assert(*it >= wp + 1);
                if (it != adjacency_lists[v].begin())
                    assert(*(it - 1) < wp + 1);
                cr_matrix(v, w) += std::distance(it, adjacency_lists[v].end());
            }
        }
    }
}

// from Dujmovic and Whitesides
void fill2(bipartite_graph &graph, matrix &cr_matrix, bool clear = true) {
    assert(graph.get_n1() == cr_matrix.get_m() &&
           cr_matrix.get_m() == cr_matrix.get_n());

    // variables
    graph.sort_adjacency_lists();
    uint64_t n0 = graph.get_n0(), n1 = graph.get_n1();
    auto adjacency_lists = graph.get_adjacency_lists();
    if (clear) memset(&cr_matrix(0, 0), 0, sizeof(uint64_t) * n1 * n1);

    // compute enhanced adjacency matrix
    // nbors(i, j) = number of nbors of i (free layer)
    // to the right of j (fixed layer)
    matrix nbors(n1, n0);
    for (uint64_t i = 0; i < n1; ++i) {
        auto it_end = adjacency_lists[i].end();
        auto it = adjacency_lists[i].begin();
        nbors(i, 0) = adjacency_lists[i].size();
        if (it != it_end && *it <= 0) {
            --nbors(i, 0);
            ++it;
        }
        for (uint64_t j = 1; j < n0; ++j) {
            nbors(i, j) = nbors(i, j - 1);
            if (it != it_end && *it <= j) {
                --nbors(i, j);
                ++it;
            }
        }
    }

    for (uint64_t v = 0; v < n1; ++v) {
        if (adjacency_lists[v].size() == 0) continue;
        uint64_t rv = adjacency_lists[v][adjacency_lists[v].size() - 1];

        for (uint64_t w = 0; w < n1; ++w) {
            if (v == w) continue;

            for (uint64_t wp : adjacency_lists[w]) {
                if (wp >= rv) break;

                cr_matrix(v, w) += nbors(v, wp);
            }
        }
    }
}

void fill_naivly(bipartite_graph &graph, matrix &matrix) {
    assert(graph.get_n1() == matrix.get_m());
    uint64_t n1 = matrix.get_m();
    auto adjacency_lists = graph.get_adjacency_lists();
    // not cache friendly
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

void print_matrix(const matrix &matrix) {
    uint64_t n1 = matrix.get_m();
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