# Cookbook demo: Writing and compiling native C++ extensions for rgraphtool
# R/C++ counterpart of graph-tool's C++ extensions cookbook (k-core decomposition,
# graph iteration, property maps, and checked vs. unchecked access).
# Runs 100% offline using Rcpp::sourceCpp() via compile_cpp_extension().

library(rgraphtool)

set.seed(2026)

cat("Public C++ header directory:\n  ", rgraphtool_include_dir(), "\n", sep = "")
cat("Bundled standalone C++ examples:\n")
cat("  kcore:           ", rgraphtool_extension_example("kcore"), "\n", sep = "")
cat("  weighted_degree: ", rgraphtool_extension_example("weighted_degree"), "\n\n", sep = "")

# 1. Build a 2,000-vertex weighted graph (core-periphery + modular structure)
n <- 2000L
m_edges <- 10000L
u <- sample.int(n, m_edges, replace = TRUE)
v <- sample.int(n, m_edges, replace = TRUE)
keep <- u != v
edges <- cbind(u[keep], v[keep])
weights <- stats::rexp(nrow(edges), rate = 1.0)

g <- graph_create(n, directed = FALSE)
graph_add_edges(g, edges, properties = list(weight = weights))
cat(sprintf("Test graph: N = %d vertices, E = %d edges\n",
            graph_num_vertices(g), graph_num_edges(g)))

# 2. Live-compile the standalone k-core and weighted-strength C++ extensions
ext_env <- new.env(parent = globalenv())
cat("Compiling kcore_extension.cpp via compile_cpp_extension()...\n")
compile_cpp_extension(file = rgraphtool_extension_example("kcore"), env = ext_env)
cat("Compiling weighted_degree_extension.cpp via compile_cpp_extension()...\n")
compile_cpp_extension(file = rgraphtool_extension_example("weighted_degree"), env = ext_env)

# 3. Verify bit-for-bit agreement with built-in functions
kc_ext_u <- ext_env$ext_kcore_decomposition(g, checked = FALSE)
kc_ext_c <- ext_env$ext_kcore_decomposition(g, checked = TRUE)
kc_ref   <- kcore_decomposition(g, checked = FALSE)
stopifnot(identical(kc_ext_u, kc_ref), identical(kc_ext_c, kc_ref))

ws_ext_u <- ext_env$ext_weighted_strength(g, "weight", checked = FALSE)
ws_ext_c <- ext_env$ext_weighted_strength(g, "weight", checked = TRUE)
ws_ref   <- vertex_strength(g, "weight", checked = FALSE)
stopifnot(max(abs(ws_ext_u - ws_ref)) < 1e-12, max(abs(ws_ext_c - ws_ref)) < 1e-12)

cat(sprintf("Max k-core shell: %d | Mean vertex strength: %.3f\n",
            max(kc_ext_u), mean(ws_ext_u)))

# 4. Compare runtime of checked vs. unchecked property maps (50 repetitions)
reps <- 50L
t_checked <- system.time({
  for (i in seq_len(reps)) ext_env$ext_kcore_decomposition(g, checked = TRUE)
})[3L]
t_unchecked <- system.time({
  for (i in seq_len(reps)) ext_env$ext_kcore_decomposition(g, checked = FALSE)
})[3L]
cat(sprintf("50 x k-core on N=%d, E=%d: checked = %.4f s | unchecked = %.4f s\n",
            n, graph_num_edges(g), t_checked, t_unchecked))

# 5. Inline custom C++ extension demonstrating triangle-wedge / two-hop edge count
inline_cpp <- '
#include <rgraphtool.h>

// [[Rcpp::export]]
Rcpp::NumericVector ext_avg_neighbor_weight(SEXP graph, std::string prop = "weight") {
    const rgraphtool::Graph& g = rgraphtool::get_graph(graph);
    const size_t n = rgraphtool::num_vertices(g);
    Rcpp::NumericVector out(static_cast<R_xlen_t>(n), 0.0);
    auto out_u = rgraphtool::vertex_property_view(out, g).get_unchecked();
    auto w_u = rgraphtool::edge_property_view(g, prop).get_unchecked();

    for (auto v : rgraphtool::vertices_range(g)) {
        double sum_w = 0.0;
        size_t cnt = 0;
        for (const auto& e : rgraphtool::out_edges_range(v, g)) {
            sum_w += w_u[e];
            ++cnt;
        }
        out_u[v] = (cnt > 0) ? (sum_w / static_cast<double>(cnt)) : 0.0;
    }
    return out;
}
'
cat("Compiling inline custom C++ extension...\n")
compile_cpp_extension(code = inline_cpp, env = ext_env)
avg_w <- ext_env$ext_avg_neighbor_weight(g, "weight")
cat(sprintf("Mean per-vertex average incident edge weight: %.4f (expected ~1.0)\n",
            mean(avg_w)))
