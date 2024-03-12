#ifndef PACE2024_COMPONENTS_HPP
#define PACE2024_COMPONENTS_HPP

#include <cassert>
#include <queue>

namespace pace2024 {

/**
 * @brief class that, given an object of `bipartite_graph`, identifies its components
 *
 * @tparam T vertex type
 */
template <typename T>
class components {
    /// @brief bipartite input graph
    const bipartite_graph<T> &graph;

    /**
     * @brief stores the component of some vertex
     * - at index 0, ..., graph.get_n_fixed() - 1 are the components of fixed vertices
     * - at index graph.get_n_fixed(), ..., graph.get_n() - 1 are the components of free vertices
     */
    std::vector<T> component_of;

    /// @brief number of components
    std::size_t nof_components{0};

   public:
    /**
     * @brief computes components of `graph` after initialization
     */
    components(const bipartite_graph<T> &graph)
        : graph(graph),
          component_of(graph.get_n(), 0) {
        identify_components();
    }

    /// @brief returns the number of components of `graph`
    std::size_t get_nof_components() {
        return nof_components;
    }

   private:
    /**
     * @brief driver method for identifying all components of `graph`
     */
    inline void identify_components() {
        // identify components of fixed layer vertices first
        const std::size_t n_fixed = graph.get_n_fixed();
        for (T v = 0; v < n_fixed; ++v) {
            if (component_of[v] == 0) {
                // v was not yet visited
                ++nof_components;
                component_of[v] = nof_components;
                identify_components_dfs(v, false);
            }
        }

        // the remaining vertices are isolated
        const std::size_t n = graph.get_n();
        for (T v = n_fixed; v < n; ++v) {
            if (component_of[v] == 0) {
                ++nof_components;
                component_of[v] = nof_components;
            }
        }
    }

    /**
     * @brief depth-first search through a component of `graph` for identifying components
     *
     * @param v current vertex
     * @param in_free_layer true if v is in free layer and false otherwise
     */
    void identify_components_dfs(const T v, const bool in_free_layer) {
        const std::vector<T> &neighbors =
            in_free_layer ? graph.get_neighbors(v) : graph.get_neighbors_fixed(v);
        for (const T &u : neighbors) {
            const T u_denormalized = in_free_layer ? u : u + graph.get_n_fixed();
            if (component_of[u_denormalized] == 0) {
                // u was not yet visited
                component_of[u_denormalized] = nof_components;
                identify_components_dfs(u, !in_free_layer);
            }
        }
    }
};

};  // namespace pace2024

#endif
