#ifndef PACE2024_MATRIX_HPP
#define PACE2024_MATRIX_HPP

#include <cstdint>
#include <type_traits>

#include "bipartite_graph.hpp"
#include "printf.hpp"

namespace pace2024 {

//
// forward declarations
//

template <typename R>
class matrix;

template <typename R>
class folded_matrix;

template <typename T, class MATRIX>
void compute_crossing_numbers_binary_search(
    const general_bipartite_graph<T> &graph, MATRIX &cr_matrix);

template <typename T, class MATRIX>
void compute_crossing_numbers_augmented_adjacency(
    const general_bipartite_graph<T> &graph, MATRIX &cr_matrix);

template <typename T, typename R>
void compute_crossing_numbers_augmented_adjacency_half(
    const general_bipartite_graph<T> &graph, folded_matrix<R> &cr_matrix);

//
// matrix types
//

/**
 * @brief generic matrix type. data is stored row-major.
 *
 * @tparam R element type
 */
template <typename R>
class matrix {
   private:
    /**
     * @brief number of rows
     */
    const std::size_t m;

    /**
     * @brief number of columns
     */
    const std::size_t n;  // number of rows and columns, respectively

    /**
     * @brief constant pointer to (mutable) data
     * matrix is stored row-major
     */
    R *const data;

   public:
    using datatype = R;

    // delete copy constructor and = operator
    matrix(const matrix &rhs) = delete;
    matrix &operator=(const matrix &rhs) = delete;
    matrix &operator=(matrix &&rhs) = delete;

    /**
     * @brief constructs the crossing number matrix from a bipartite graph
     * that resembles a one-sided crossing minimization instance.
     *
     * @tparam T vertex type
     * @param graph
     */
    template <typename T>
    matrix(const general_bipartite_graph<T> &graph)
        : m(graph.get_n1()),  //
          n(graph.get_n1()),
          data(new R[m * n]()) {
        // the additional () above initializes memory to 0
        assert(graph.get_n1() > 0);
        compute_crossing_numbers_binary_search(graph, *this);
    }

    /**
     * @brief initializes a (m x n) matrix with 0s
     *
     * @param m
     * @param n
     */
    matrix(const std::size_t m, const std::size_t n)
        : m(m), n(n), data(new R[m * n]()) {
        // the additional () initializes memory to 0
        assert(m > 0);
        assert(n > 0);
    }

    /**
     * @brief deletes memory behind data pointer
     */
    ~matrix() { delete[] data; }

    /**
     * @brief overrides the matrix with 0s
     */
    void clear() { memset(data, 0, sizeof(R) * m * n); }

    /**
     * @brief returns the number of rows
     *
     * @return std::size_t
     */
    std::size_t get_m() const { return m; }

    /**
     * @brief returns the number of columns
     *
     * @return std::size_t
     */
    std::size_t get_n() const { return n; }

    /**
     * @brief returns reference to element at row i and column j
     *
     * @param i
     * @param j
     * @return R&
     */
    R &operator()(std::size_t i, std::size_t j) { return data[i * n + j]; }

    /**
     * @brief returns constant reference to element at row i and column j
     *
     * @param i
     * @param j
     * @return const R&
     */
    const R &operator()(std::size_t i, std::size_t j) const {
        return data[i * n + j];
    }
};

/**
 * @brief special matrix for storing crossing number matrix of an instance
 * of one-sided crossing minimization.
 * data layout is chosen s.t. iterating it is cache-friendly.
 * essentially it is a square matrix w/o diagonal elements.
 *
 * @tparam R element type
 */
template <typename R>
class folded_matrix {
   private:
    /**
     * @brief number of vertices in the free partite set
     */
    const std::size_t n1;

    /**
     * @brief n2 = n1 * (n1 - 1)
     */
    const std::size_t n2;

    /**
     * @brief constant pointer to (mutable) data
     * storage order: c_01, c_10, c_02, c_20, ..., c_12, c_21, ...
     * no diagonal elements
     */
    R *const data;

   public:
    using datatype = R;

    // delete copy constructor and = operator
    folded_matrix(const folded_matrix &rhs) = delete;
    folded_matrix &operator=(const folded_matrix &rhs) = delete;
    folded_matrix &operator=(folded_matrix &&rhs) = delete;

    /**
     * @brief constructs the crossing number matrix from a bipartite graph
     * that resembles a one-sided crossing minimization instance.
     *
     * @tparam T vertex type
     * @param graph
     */
    template <typename T>
    folded_matrix(const general_bipartite_graph<T> &graph)
        : n1(graph.get_n1()), n2(n1 * (n1 - 1)), data(new R[n2]()) {
        // the additional () above initializes memory to 0
        assert(n1 > 1);
        compute_crossing_numbers_augmented_adjacency_half(graph, *this);
    }

    /**
     * @brief constructs folded_matrix initialized with 0s
     *
     * @param n1
     */
    folded_matrix(const std::size_t n1)
        : n1(n1), n2(n1 * (n1 - 1)), data(new R[n2]()) {
        // the additional () above initializes memory to 0
        assert(n1 > 1);
    }

    /**
     * @brief deletes memory behind pointer data
     */
    ~folded_matrix() { delete[] data; }

    /**
     * @brief overwrites data with 0s
     */
    void clear() { memset(data, 0, sizeof(R) * n2); }

    /**
     * @brief returns number of rows (matrix is square)
     *
     * @return std::size_t
     */
    std::size_t get_m() const { return n1; }

    /**
     * @brief returns number of columns (matrix is square)
     *
     * @return std::size_t
     */
    std::size_t get_n() const { return n1; }

    /**
     * @brief converts indices i (row) and j (column) to the
     * correct position of the corresponding element in data.
     *
     * @param i
     * @param j
     * @return std::size_t
     */
    inline std::size_t get_index(const std::size_t i,
                                 const std::size_t j) const {
        assert(i < n1);
        assert(j < n1);
        assert(i != j);
        std::size_t index;
        if (i < j) {
            std::size_t offset = n2 - (n1 - i) * (n1 - i - 1);
            index = offset + 2 * (j - i - 1);
        } else {
            std::size_t offset = n2 - (n1 - j) * (n1 - j - 1);
            index = offset + 2 * (i - j - 1) + 1;
        }
        assert(index < n2);
        return index;
    }

    /**
     * @brief returns reference to the element in row i and column j.
     * we require i != j.
     *
     * @param i
     * @param j
     * @return T&
     */
    R &operator()(const std::size_t i, const std::size_t j) {
        return data[get_index(i, j)];
    }

    /**
     * @brief returns constant reference to the element in row i and column j.
     * we require i != j.
     *
     * @param i
     * @param j
     * @return const T&
     */
    const R &operator()(const std::size_t i, const std::size_t j) const {
        return data[get_index(i, j)];
    }
};

using uint64_folded_matrix = folded_matrix<uint64_t>;
using uint32_folded_matrix = folded_matrix<uint32_t>;
using uint16_folded_matrix = folded_matrix<uint16_t>;

//
// compute crossing numbers
//

/**
 * @brief fills the matrix with crossing numbers.
 * adapted from Dujmovic's and Whitesides' paper https://doi.org/10.1007/s00453-004-1093-2
 *
 * @tparam T vertex type
 * @tparam MATRIX a matrix type
 * @param graph
 * @param cr_matrix
 */
template <typename T, class MATRIX>
void compute_crossing_numbers_binary_search(
    const general_bipartite_graph<T> &graph, MATRIX &cr_matrix) {
    assert(graph.get_n1() == cr_matrix.get_m() && cr_matrix.get_m() == cr_matrix.get_n());
    std::size_t n1 = graph.get_n1();
    auto &adjacency_lists = graph.get_adjacency_lists();

    for (T v = 0; v < n1; ++v) {
        if (adjacency_lists[v].size() == 0) continue;
        T rv = adjacency_lists[v][adjacency_lists[v].size() - 1];

        for (T w = 0; w < n1; ++w) {
            if (v == w) continue;

            cr_matrix(v, w) = 0;
            for (T wp : adjacency_lists[w]) {
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

/**
 * @brief fills the matrix with the crossing numbers.
 * from Dujmovic and Whitesides https://doi.org/10.1007/s00453-004-1093-2
 *
 * @tparam T vertex type
 * @tparam MATRIX a matrix type
 * @param graph
 * @param cr_matrix
 */
template <typename T, class MATRIX>
void compute_crossing_numbers_augmented_adjacency(
    const general_bipartite_graph<T> &graph, MATRIX &cr_matrix) {
    assert(graph.get_n1() == cr_matrix.get_m() && cr_matrix.get_m() == cr_matrix.get_n());

    // variables
    std::size_t n0 = graph.get_n0(), n1 = graph.get_n1();
    auto &adjacency_lists = graph.get_adjacency_lists();

    // compute enhanced adjacency matrix
    // nbors(i, j) = number of nbors of i (free layer)
    // to the right of j (fixed layer)
    matrix<typename MATRIX::datatype> nbors(n1, n0);
    for (T i = 0; i < n1; ++i) {
        auto it_end = adjacency_lists[i].end();
        auto it = adjacency_lists[i].begin();
        nbors(i, 0) = adjacency_lists[i].size();
        if (it != it_end && *it <= 0) {
            --nbors(i, 0);
            ++it;
        }
        for (T j = 1; j < n0; ++j) {
            nbors(i, j) = nbors(i, j - 1);
            if (it != it_end && *it <= j) {
                --nbors(i, j);
                ++it;
            }
        }
    }

    for (T v = 0; v < n1; ++v) {
        if (adjacency_lists[v].size() == 0) continue;
        T rv = adjacency_lists[v][adjacency_lists[v].size() - 1];

        for (T w = 0; w < n1; ++w) {
            if (v == w) continue;

            cr_matrix(v, w) = 0;
            for (T wp : adjacency_lists[w]) {
                if (wp >= rv) break;

                cr_matrix(v, w) += nbors(v, wp);
            }
        }
    }
}

/**
 * @brief fills a folded_matrix with crossing numbers.
 * adapted from Dujmovic and Whitesides https://doi.org/10.1007/s00453-004-1093-2
 *
 * @tparam T vertex type
 * @param graph
 * @param cr_matrix
 */
template <typename T, typename R>
void compute_crossing_numbers_augmented_adjacency_half(
    const general_bipartite_graph<T> &graph, folded_matrix<R> &cr_matrix) {
    assert(graph.get_n1() == cr_matrix.get_n());

    // variables
    std::size_t n0 = graph.get_n0(), n1 = graph.get_n1();
    auto &adjacency_lists = graph.get_adjacency_lists();

    // compute enhanced adjacency matrix
    // nbors(i, j) = number of nbors of i (free layer)
    // to the right of j (fixed layer)
    matrix<R> nbors(n1, n0);
    for (T i = 0; i < n1; ++i) {
        auto it_end = adjacency_lists[i].end();
        auto it = adjacency_lists[i].begin();
        nbors(i, 0) = adjacency_lists[i].size();
        if (it != it_end && *it <= 0) {
            --nbors(i, 0);
            ++it;
        }
        for (T j = 1; j < n0; ++j) {
            nbors(i, j) = nbors(i, j - 1);
            if (it != it_end && *it <= j) {
                --nbors(i, j);
                ++it;
            }
        }
    }

    for (T v = 0; v < n1; ++v) {
        const std::size_t deg_v = adjacency_lists[v].size();
        if (deg_v == 0) continue;
        const T rv = adjacency_lists[v][deg_v - 1];

        for (T w = v + 1; w < n1; ++w) {
            cr_matrix(v, w) = 0;

            std::size_t deg_w = adjacency_lists[w].size();
            if (deg_w == 0) continue;

            for (T wp : adjacency_lists[w]) {
                if (wp >= rv) break;

                cr_matrix(v, w) += nbors(v, wp);
            }

            R nof_common_nbors = 0;
            for (T i = 0, j = 0; i < deg_v && j < deg_w;) {
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

namespace test {

/**
 * @brief fills the matrix with all crossing numbers by
 * considering the ends of all pairs of edges.
 *
 * @tparam T vertex type
 * @tparam R element type
 * @param graph
 * @param matrix
 */
template <typename T, typename R>
void compute_crossing_numbers_naivly(const general_bipartite_graph<T> &graph,
                                     folded_matrix<R> &matrix) {
    assert(graph.get_n1() == matrix.get_m());

    const std::size_t n1 = matrix.get_m();
    auto &adjacency_lists = graph.get_adjacency_lists();
    for (T i = 0; i < n1 - 1; ++i) {
        for (T j = i + 1; j < n1; ++j) {
            for (const T &x : adjacency_lists[i]) {
                for (const T &y : adjacency_lists[j]) {
                    if (x < y)
                        ++matrix(j, i);
                    else if (x > y)
                        ++matrix(i, j);
                }
            }
        }
    }
}

/**
 * @brief prints the matrix.
 * not cache-friendly.
 *
 * @tparam R element type
 * @param matrix
 */
template <typename R>
void print_matrix(const matrix<R> &matrix) {
    const std::size_t m = matrix.get_m();
    const std::size_t n = matrix.get_n();
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            // Bossert's printf does not care for the right wildcard <3
            fmt::printf("%5s", matrix(i, j));
        }
        fmt::printf("\n");
    }
}

};  // namespace test

}  // namespace pace2024

#endif