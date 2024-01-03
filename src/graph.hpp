#ifndef PACE2024_GRAPH_HPP
#define PACE2024_GRAPH_HPP

#include <cstdint>
#include <vector>

namespace pace2024 {

template <typename T>
class general_graph {
   private:
    std::vector<std::vector<T>> adjacency_lists;

   public:
    using datatype = T;

    general_graph(const general_graph &rhs) = delete;
    general_graph &operator=(const general_graph &rhs) = delete;
    general_graph &operator=(general_graph &&rhs) = delete;

    general_graph(const std::size_t n) : adjacency_lists(n) {}

    /**
     * @brief add arc (u,v)
     * 
     * @param u 
     * @param v 
     */
    void add_edge(const T u, const T v) {
        adjacency_lists[u].emplace_back(v);
    }

    std::size_t get_n() const {
        return adjacency_lists.size();
    }

    const std::vector<T>& get_adjacency_list(const T v) const {
        return adjacency_lists[v];
    }
};

};  // namespace pace2024

#endif