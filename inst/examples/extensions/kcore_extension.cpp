// Standalone C++ extension example for rgraphtool: k-core decomposition
// R/C++ counterpart of graph-tool's k-core C++ extension tutorial.
//
// Compile from R with:
//   rgraphtool::compile_cpp_extension(rgraphtool::rgraphtool_extension_example("kcore"))
// or (when rgraphtool is installed):
//   Rcpp::sourceCpp(rgraphtool::rgraphtool_extension_example("kcore"))

// [[Rcpp::depends(rgraphtool, BH)]]
// [[Rcpp::plugins(cpp23)]]

#include <rgraphtool.h>
#include <vector>

// [[Rcpp::export]]
Rcpp::IntegerVector ext_kcore_decomposition(SEXP graph, bool checked = false) {
    const rgraphtool::Graph& g = rgraphtool::get_graph(graph);
    const size_t n = rgraphtool::num_vertices(g);
    Rcpp::IntegerVector core_vec(static_cast<R_xlen_t>(n));
    if (n == 0) return core_vec;

    // Property maps: checked by default, convertible to zero-overhead unchecked views
    auto core_map = rgraphtool::vertex_property_view(core_vec, g);
    auto deg_map = rgraphtool::make_vertex_property<size_t>(g, 0);
    auto pos_map = rgraphtool::make_vertex_property<size_t>(g, 0);

    auto run_kcore = [&](auto& core, auto& deg, auto& pos) {
        std::vector<std::vector<rgraphtool::Vertex>> bins;
        for (auto v : rgraphtool::vertices_range(g)) {
            size_t k = rgraphtool::total_degree(v, g);
            deg[v] = k;
            if (k >= bins.size()) bins.resize(k + 1);
            pos[v] = bins[k].size();
            bins[k].push_back(v);
        }
        for (size_t k = 0; k < bins.size(); ++k) {
            while (!bins[k].empty()) {
                auto v = bins[k].back();
                bins[k].pop_back();
                core[v] = static_cast<int>(k);
                for (auto u : rgraphtool::all_neighbors_range(v, g)) {
                    if (deg[u] > deg[v]) {
                        size_t ku = deg[u];
                        auto w = bins[ku].back();
                        size_t pu = pos[u];
                        bins[ku][pu] = w;
                        pos[w] = pu;
                        bins[ku].pop_back();
                        pos[u] = bins[ku - 1].size();
                        bins[ku - 1].push_back(u);
                        deg[u] = ku - 1;
                    }
                }
            }
        }
    };

    if (checked) {
        run_kcore(core_map, deg_map, pos_map);
    } else {
        auto core_u = core_map.get_unchecked();
        auto deg_u = deg_map.get_unchecked();
        auto pos_u = pos_map.get_unchecked();
        run_kcore(core_u, deg_u, pos_u);
    }
    return core_vec;
}
