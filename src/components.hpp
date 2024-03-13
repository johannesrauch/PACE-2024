#ifndef PACE2024_COMPONENTS_HPP
#define PACE2024_COMPONENTS_HPP

#include <cassert>
#include <cstdint>
#include <queue>
#include <unordered_map>

#include "bipartite_graph.hpp"

namespace pace2024 {

/**
 * @brief class that, given an object of `bipartite_graph`, identifies its components
 *
 * @tparam T vertex type
 */
template <typename T>
class components {
    /// @brief bipartite input graph
    const bipartite_graph<T>& graph;

    /**
     * @brief stores the component of some vertex
     * - at index 0, ..., graph.get_n_fixed() - 1 are the components of fixed vertices
     * - at index graph.get_n_fixed(), ..., graph.get_n() - 1 are the components of free vertices
     */
    std::vector<T> component_of;

    /// @brief number of components
    std::size_t nof_components{0};

    /// @brief stores the components of `graph`
    std::vector<bipartite_graph<T>> components_;

    std::vector<std::vector<T>> original_free_vertices;

    bool built{false};

   public:
    /**
     * @brief computes components of `graph` after initialization
     */
    components(const bipartite_graph<T>& graph)
        : graph(graph),
          component_of(graph.get_n(), 0) {
        identify();
    }

    //
    // for each support
    //

    void build() {
        if (built) return;
        built = true;

        // add fixed layer vertices to components
        components_.reserve(nof_components);
        for (std::size_t i = 0; i < nof_components; ++i) components_.emplace_back();  // resize does not call standard ctor
        std::vector<std::unordered_map<T, T>> maps(nof_components);
        for (T v = 0; v < graph.get_n_fixed(); ++v) {
            assert(component_of[v] >= 1);
            const T c = component_of[v] - 1;
            maps[c][v] = components_[c].add_fixed_vertex();
        }

        // add free layer vertices to components
        original_free_vertices.resize(nof_components);
        for (T v = graph.get_n_fixed(); v < graph.get_n(); ++v) {
            assert(component_of[v] >= 1);
            const T c = component_of[v] - 1;
            maps[c][v] = components_[c].add_free_vertex();
            original_free_vertices[c].emplace_back(v - graph.get_n_fixed());
        }

        // add edges
        for (std::size_t c = 0; c < nof_components; ++c) {
            assert(!maps[c].empty());
            const auto [v_index, v_component] = *maps[c].begin();
            (void)v_component;
            build_bfs(v_index, maps[c], components_[c]);
        }
    }

    //
    // getter
    //

    const bipartite_graph<T>& get_component(const std::size_t c) const {
        assert(c < nof_components);
        if (nof_components == 1) {
            return graph;
        } else {
            assert(built);
            return components_[c];
        }
    }

    /// @brief returns the number of components of `graph`
    std::size_t get_nof_components() {
        return nof_components;
    }

   private:
    inline void build_bfs(const T v_index, std::unordered_map<T, T> map, bipartite_graph<T>& component) {
        std::queue<std::pair<T, bool>> q;
        const bool v_in_free_layer = v_index >= graph.get_n_fixed();
        const T v = v_index - (v_in_free_layer ? graph.get_n_fixed() : 0);
        q.emplace(v, v_in_free_layer);
        std::vector<uint8_t> visited(graph.get_n(), 0);

        while (!q.empty()) {
            const auto [w, in_free_layer] = q.front();
            q.pop();
            const T w_index = w + (in_free_layer ? graph.get_n_fixed() : 0);
            const auto& neighbors = in_free_layer ? graph.get_neighbors_of_free(w) : graph.get_neighbors_of_fixed(w);

            for (const T& u : neighbors) {
                const T u_index = u + (in_free_layer ? 0 : graph.get_n_fixed());
                if (visited[u_index] == 0) {
                    visited[u_index] = 1;
                    q.emplace(u, !in_free_layer);
                }
                if (visited[u_index] == 1) {
                    if (in_free_layer) {
                        component.add_edge(map[u_index], map[w_index]);
                    } else {
                        component.add_edge(map[w_index], map[u_index]);
                    }
                }
            }
            visited[w_index] = 2;
        }
    }

    /**
     * @brief driver method for identifying and building all components of `graph`
     */
    inline void identify() {
        // identify components of fixed layer vertices first
        const std::size_t n_fixed = graph.get_n_fixed();
        for (T v = 0; v < n_fixed; ++v) {
            if (component_of[v] == 0) {
                // v was not yet visited
                ++nof_components;
                component_of[v] = nof_components;
                identify_bfs(v);
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
     * @brief breadth-first search through a component of `graph` to identify its vertices
     *
     * @param v some vertex of the fixed layer
     */
    void identify_bfs(const T v_fixed) {
        std::queue<std::pair<T, bool>> q;
        q.emplace(v_fixed, false);

        while (!q.empty()) {
            const auto [w, in_free_layer] = q.front();
            q.pop();
            const auto& neighbors = in_free_layer ? graph.get_neighbors_of_free(w) : graph.get_neighbors_of_fixed(w);

            for (const T& u : neighbors) {
                const T u_index = u + (in_free_layer ? 0 : graph.get_n_fixed());
                const bool not_yet_visited = component_of[u_index] == 0;
                if (not_yet_visited) {
                    component_of[u_index] = nof_components;
                    q.emplace(u, !in_free_layer);
                }
            }
        }
    }
};

};  // namespace pace2024

#endif
