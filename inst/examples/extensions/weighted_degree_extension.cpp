// Standalone C++ extension example for rgraphtool: weighted vertex strength
// Demonstrates edge property map access, vertex property map writing, and
// checked vs. unchecked property map views via <rgraphtool.h>.
//
// Compile from R with:
//   rgraphtool::compile_cpp_extension(
//     rgraphtool::rgraphtool_extension_example("weighted_degree")
//   )

// [[Rcpp::depends(rgraphtool, BH)]]
// [[Rcpp::plugins(cpp23)]]

#include <rgraphtool.h>
#include <string>

// [[Rcpp::export]]
Rcpp::NumericVector ext_weighted_strength(SEXP graph,
                                          std::string edge_prop = "weight",
                                          bool checked = false) {
    const rgraphtool::Graph& g = rgraphtool::get_graph(graph);
    const size_t n = rgraphtool::num_vertices(g);
    Rcpp::NumericVector strength_vec(static_cast<R_xlen_t>(n), 0.0);
    if (n == 0) return strength_vec;

    auto strength_map = rgraphtool::vertex_property_view(strength_vec, g);
    auto weight_map = rgraphtool::edge_property_view(g, edge_prop);

    auto run_strength = [&](auto& strength, const auto& weight) {
        for (auto v : rgraphtool::vertices_range(g)) {
            double total = 0.0;
            for (const auto& e : rgraphtool::all_edges_range(v, g)) {
                total += weight[e];
            }
            strength[v] = total;
        }
    };

    if (checked) {
        run_strength(strength_map, weight_map);
    } else {
        auto strength_u = strength_map.get_unchecked();
        auto weight_u = weight_map.get_unchecked();
        run_strength(strength_u, weight_u);
    }
    return strength_vec;
}
