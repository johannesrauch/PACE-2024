#include "log/debug_printf.hpp"

namespace pace {

namespace test {

void printf_info(const branch_and_cut_info& info_br,
                        const highs_lp_info& info_lp) {
    if (info_br.n_iters % 50 == 0) {
        fmt::printf(
            "%11s%11s%11s|%11s%11s%11s%11s%11s%11s|%11s%11s%11s%11s|%11s%11s%"
            "11s%11s|%11s%11s%11s%11s%11s|%11s\n",                         //
            "walltime", "lb", "ub",                                        //
            "iter", "obj val", "t simplex", "n rows", "n nodes", "depth",  //
            "n simplex", "n splx avg", "n cold",
            "was warm",                                     //
            "n bucket", "n added", "n delete", "n spared",  //
            "uvw iter", "on node", "u", "v", "w",           //
            "rh conf");
    }
    fmt::printf(
        "%11.1f%11u%11u|%11u%11.1f%11.1f%11u%11u%11u|%11u%11.1f%11u%11u|%11u%"
        "11u%11u%11u|%11u%11u%11u%11u%11u|%11.3f\n",  //
        elapsed_walltime_in_s(),                      //
        info_br.lower_bound,                          //
        info_br.upper_bound,                          //
        info_br.n_iters,                              //
        info_lp.objective_value,                      //
        info_lp.t_simplex,                            //
        info_lp.n_rows,                               //
        info_br.n_search_nodes,                       //
        info_br.depth,                                //
        info_lp.n_simplex_iters,                      //
        info_lp.n_avg_simplex_iters,                  //
        info_lp.n_simplex_coldstart_iters,            //
        info_lp.was_warmstart,                        //
        info_lp.n_bucket_entries,                     //
        info_lp.n_added_rows,                         //
        info_lp.n_deleted_rows,                       //
        info_lp.n_delete_rows_spared,                 //
        info_lp.n_3cycle_iters,                       //
        info_br.n_iter_3cycle_current_node,           //
        info_lp.u_old,                                //
        info_lp.v_old,                                //
        info_lp.w_old,                                //
        info_br.relax_h_confidence);
}

void printf_lpinfo_line() {
    fmt::printf(
        "%11s|%11s%11s%11s%11s%11s|%11s%11s%11s%11s|%11s%11s%11s%11s\n",  //
        "walltime",                                                       //
        "n iters", "obj val", "n rows", "n splx", "n splx avg",           //
        "n bucket", "n added", "n delete", "n spared",                    //
        "uvw iter", "u", "v", "w");
    std::cout << std::flush;
}

void printf_lpinfo(const highs_lp_info& info) {
    fmt::printf(
        "%11.1f|%11u%11.1f%11u%11u%11.1f|%11u%11u%11u%11u|%11u%11u%11u%"
        "11u\n",                    //
        elapsed_walltime_in_s(),    //
        info.n_solve_iters,         //
        info.objective_value,       //
        info.n_rows,                //
        info.n_simplex_iters,       //
        info.n_avg_simplex_iters,   //
        info.n_bucket_entries,      //
        info.n_added_rows,          //
        info.n_deleted_rows,        //
        info.n_delete_rows_spared,  //
        info.n_3cycle_iters,        //
        info.u_old, info.v_old, info.w_old);
    std::cout << std::flush;
}

};  // namespace test

};  // namespace pace
