#ifndef PACE2024_MATRIX_HPP
#define PACE2024_MATRIX_HPP

#include <cstdint>

#include "bipartite_graph.hpp"
#include "printf.hpp"

namespace pace2024 {

template <typename T>
class matrix;

template <typename T>
class folded_square_matrix;

template <typename T>
void fill(bipartite_graph &graph, matrix<T> &matrix);

template <class MATRIX>
void fill2(bipartite_graph &graph, MATRIX &cr_matrix);

template <typename T>
void fill(bipartite_graph &graph, folded_square_matrix<T> &cr_matrix);

template <typename T>
class matrix {
   private:
    const std::size_t m, n;  // number of rows and columns, respectively
    T *const data;           // row-major

   public:
    using datatype = T;

    matrix(const matrix &rhs) = delete;
    matrix &operator=(const matrix &rhs) = delete;
    matrix &operator=(matrix &&rhs) = delete;

    matrix(bipartite_graph &graph)
        : m(graph.get_n1()),  //
          n(graph.get_n1()),
          data(new T[m * n]()) {
        // the additional () above initializes memory to 0
        assert(graph.get_n1() > 0);
        fill2(graph, *this);
    }

    matrix(const uint64_t m, const uint64_t n)
        : m(m), n(n), data(new T[m * n]()) {
        assert(m > 0);
        assert(n > 0);
    }  // the additional () initializes memory to 0

    ~matrix() { delete[] data; }

    void clear() { memset(data, 0, sizeof(T) * m * n); }

    std::size_t get_m() const { return m; }

    std::size_t get_n() const { return n; }

    T &operator()(std::size_t i, std::size_t j) { return data[i * n + j]; }

    const T &operator()(std::size_t i, std::size_t j) const {
        return data[i * n + j];
    }
};

template <typename T>
class folded_square_matrix {
   private:
    const std::size_t n1,  // number of vertices in the free partite set
        n2;                // = n1 * (n1 - 1)
    T *const data;  // storage order: c_01, c_10, c_02, c_20, ..., c_12, c_21,
                    // ... no diagonal elements

   public:
    using datatype = T;

    folded_square_matrix(const folded_square_matrix &rhs) = delete;
    folded_square_matrix &operator=(const folded_square_matrix &rhs) = delete;
    folded_square_matrix &operator=(folded_square_matrix &&rhs) = delete;

    folded_square_matrix(bipartite_graph &graph)
        : n1(graph.get_n1()), n2(n1 * (n1 - 1)), data(new T[n2]()) {
        // the additional () above initializes memory to 0
        assert(n1 > 1);
        fill(graph, *this);
    }

    folded_square_matrix(const std::size_t n1)
        : n1(n1), n2(n1 * (n1 - 1)), data(new T[n2]()) {
        // the additional () above initializes memory to 0
        assert(n1 > 1);
    }

    ~folded_square_matrix() { delete[] data; }

    void clear() { memset(data, 0, sizeof(T) * n2); }

    std::size_t get_m() const { return n1; }

    std::size_t get_n() const { return n1; }

    T &operator()(const std::size_t i, const std::size_t j) {
        assert(i < n1 && j < n1 && i != j);
        std::size_t index;
        if (i < j) {
            std::size_t offset = n2 - (n1 - i) * (n1 - i - 1);
            index = offset + 2 * (j - i - 1);
        } else {
            std::size_t offset = n2 - (n1 - j) * (n1 - j - 1);
            index = offset + 2 * (i - j - 1) + 1;
        }
        // std::cout << "n1=" << n1 << ",n2=" << n2 << ",i=" << i << ",j=" << j
        //          << ",idx=" << index << std::endl;
        assert(index < n2);
        return data[index];
    }

    const T &operator()(const std::size_t i, const std::size_t j) const {
        assert(i < n1 && j < n1 && i != j);
        std::size_t index;
        if (i < j) {
            std::size_t offset = n2 - (n1 - i) * (n1 - i - 1);
            index = offset + 2 * (j - i - 1);
        } else {
            std::size_t offset = n2 - (n1 - j) * (n1 - j - 1);
            index = offset + 2 * (i - j - 1) + 1;
        }
        assert(index < n2);
        return data[index];
    }
};

// adapted from Dujmovic and Whitesides
template <typename T>
void fill(bipartite_graph &graph, matrix<T> &cr_matrix) {
    assert(graph.get_n1() == cr_matrix.get_m());
    graph.sort_adjacency_lists();
    uint64_t n1 = graph.get_n1();
    auto adjacency_lists = graph.get_adjacency_lists();

    for (uint64_t v = 0; v < n1; ++v) {
        if (adjacency_lists[v].size() == 0) continue;
        uint64_t rv = adjacency_lists[v][adjacency_lists[v].size() - 1];

        for (uint64_t w = 0; w < n1; ++w) {
            if (v == w) continue;

            cr_matrix(v, w) = 0;
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
template <class MATRIX>
void fill2(bipartite_graph &graph, MATRIX &cr_matrix) {
    assert(graph.get_n1() == cr_matrix.get_m() &&
           cr_matrix.get_m() == cr_matrix.get_n());

    // variables
    graph.sort_adjacency_lists();
    uint64_t n0 = graph.get_n0(), n1 = graph.get_n1();
    auto adjacency_lists = graph.get_adjacency_lists();

    // compute enhanced adjacency matrix
    // nbors(i, j) = number of nbors of i (free layer)
    // to the right of j (fixed layer)
    matrix<typename MATRIX::datatype> nbors(n1, n0);
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

            cr_matrix(v, w) = 0;
            for (uint64_t wp : adjacency_lists[w]) {
                if (wp >= rv) break;

                cr_matrix(v, w) += nbors(v, wp);
            }
        }
    }
}

// adapted from Dujmovic and Whitesides
template <typename T>
void fill(bipartite_graph &graph, folded_square_matrix<T> &cr_matrix) {
    assert(graph.get_n1() == cr_matrix.get_n());

    // variables
    graph.sort_adjacency_lists();
    uint64_t n0 = graph.get_n0(), n1 = graph.get_n1();
    auto adjacency_lists = graph.get_adjacency_lists();

    // compute enhanced adjacency matrix
    // nbors(i, j) = number of nbors of i (free layer)
    // to the right of j (fixed layer)
    matrix<T> nbors(n1, n0);
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
        const std::size_t deg_v = adjacency_lists[v].size();
        if (deg_v == 0) continue;
        const uint64_t rv = adjacency_lists[v][deg_v - 1];

        for (uint64_t w = v + 1; w < n1; ++w) {
            cr_matrix(v, w) = 0;

            std::size_t deg_w = adjacency_lists[w].size();
            if (deg_w == 0) continue;

            for (uint64_t wp : adjacency_lists[w]) {
                if (wp >= rv) break;

                cr_matrix(v, w) += nbors(v, wp);
            }

            T nof_common_nbors = 0;
            for (std::size_t i = 0, j = 0; i < deg_v && j < deg_w;) {
                if (adjacency_lists[v][i] == adjacency_lists[w][j]) {
                    ++nof_common_nbors;
                    ++i;
                    ++j;
                } else if (adjacency_lists[v][i] < adjacency_lists[w][j]) {
                    ++i;
                } else {
                    ++j;
                }
            }
            cr_matrix(w, v) =
                deg_v * deg_w - nof_common_nbors - cr_matrix(v, w);
        }
    }
}

template <typename T>
void fill_naivly(bipartite_graph &graph, folded_square_matrix<T> &matrix) {
    assert(graph.get_n1() == matrix.get_m());
    uint64_t n1 = matrix.get_m();
    auto adjacency_lists = graph.get_adjacency_lists();
    // not cache friendly for matrix<T>
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

template <typename T>
void print_matrix(const matrix<T> &matrix) {
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