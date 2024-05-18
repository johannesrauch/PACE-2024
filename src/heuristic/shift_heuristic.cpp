#include "heuristic/shift_heuristic.hpp"

namespace pace {

crossing_number_t shift_heuristic::operator()(std::vector<vertex_t> &ordering,
                                              crossing_number_t n_crossings) {
    assert(ordering.size() == n_free);
    PACE_DEBUG_PRINTF("start shift heuristic, shift_length=%11u\n", params.shift_length);
    if (lower_bound() >= n_crossings) return n_crossings;

    assert(ordering.size() == graph.get_n_free());
    std::copy(ordering.begin(), ordering.end(), ordering_.begin());
    uint32_t go_on = 1;

    for (std::size_t i = 0;  //
         i < params.n_iterations && lower_bound() < n_crossings && go_on &&
         !tle;
         ++i) {
        go_on = false;
        for (typename std::list<vertex_t>::iterator it = ordering_.begin();
             it != ordering_.end() && !tle;) {
            const auto [go_on_, it_] = improve(it, n_crossings);
            go_on += go_on_;
            it = it_;
        }
        PACE_DEBUG_PRINTF("%11s=%11s\n", "n improve", go_on);
    }

    std::copy(ordering_.begin(), ordering_.end(), ordering.begin());
    assert(ordering.size() == graph.get_n_free());
    assert(number_of_crossings(graph, ordering) == n_crossings);
    update_ordering(ordering, n_crossings);

    PACE_DEBUG_PRINTF("end   shift heuristic\n");
    return n_crossings;
}

std::pair<bool, typename std::list<vertex_t>::iterator>  //
shift_heuristic::improve(                                //
    typename std::list<vertex_t>::iterator i, uint32_t &n_crossings) {
    // try left shifts
    uint32_t c_old{0}, c_new{0}, improvement = 0;
    typename std::list<vertex_t>::iterator k = i, l = i;
    for (std::size_t j = 1; j <= params.shift_length; ++j) {
        if (k == ordering_.begin()) break;
        --k;

        const auto [c1, c2] = cr_numbers(*k, *i);
        c_old += c1;
        c_new += c2;

        if (c_new < c_old && c_old - c_new > improvement) {
            improvement = c_old - c_new;
            l = k;
        }
    }

    // try right shifts
    c_old = 0;
    c_new = 0;
    k = i;
    bool right_shift = false;
    for (std::size_t j = 1; j <= params.shift_length; ++j) {
        ++k;
        if (k == ordering_.end()) break;

        const auto [c1, c2] = cr_numbers(*i, *k);
        c_old += c1;
        c_new += c2;

        if (c_new < c_old && c_old - c_new > improvement) {
            improvement = c_old - c_new;
            l = k;
            right_shift = true;
        }
    }

    if (improvement == 0) return std::make_pair(false, ++i);

    // move element at i in front of element at l
    if (right_shift) ++l;
    ordering_.insert(l, *i);
    typename std::list<vertex_t>::iterator it = ordering_.erase(i);

    assert(ordering_.size() == graph.get_n_free());
    assert(n_crossings >= improvement);
    n_crossings -= improvement;
    assert(number_of_crossings(graph, std::vector<vertex_t>{ordering_.begin(),
                                                            ordering_.end()}) ==
           n_crossings);

    return std::make_pair(true, it);
}

}  // namespace pace
