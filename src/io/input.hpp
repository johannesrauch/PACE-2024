#ifndef PACE_IO_INPUT_HPP
#define PACE_IO_INPUT_HPP

#include <limits>
#include <map>
#include <vector>

#include "io/parse_input.hpp"
#include "log/debug_printf.hpp"
#include "model/bipartite_graph.hpp"

namespace pace {

/**
 * @brief general input that is able to split into (smaller) instances
 */
class input {
    bipartite_graph graph;
    std::vector<std::unique_ptr<bipartite_graph>> subgraph_ptrs;

    std::vector<vertex_t> free_layer;
    std::vector<vertex_t> borders;

    bool first_graph_empty{false};
    bool tried_split{false};

   public:
    const fs::path filepath;

    input();

    input(const fs::path filepath);

    // delete copy constructor and assignment
    input(const input &other) = delete;
    input &operator=(const input &other) = delete;

    const bipartite_graph &get_graph() const { return graph; }

    std::size_t get_n_subgraphs() {
        if (!tried_split) try_split();
        return subgraph_ptrs.size();
    }

    const bipartite_graph &get_subgraph(const std::size_t i) {
        if (!tried_split) try_split();
        assert(i < get_n_subgraphs());
        return *subgraph_ptrs[i];
    }

    bool is_first_graph_empty() {
        if (!tried_split) try_split();
        return first_graph_empty;
    }

    void lift_ordering(const std::size_t i,
                       const std::vector<vertex_t> &subordering,
                       std::vector<vertex_t> &ordering);

    bool try_split();

    bool compare(const vertex_t &u, const vertex_t &v) const;

   private:
    inline void add_subgraph(const std::size_t begin, const std::size_t end);

    inline void sort_free_layer();
};

};  // namespace pace

#endif
