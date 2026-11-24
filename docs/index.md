# rgraphtool

`rgraphtool` is a standalone **native R/C++ port** of the Python
[graph-tool](https://graph-tool.skewed.de/) network analysis library
(originally developed by Tiago P. Peixoto). It embeds `graph-tool`’s C++
adjacency data structure (`boost::adj_list`) and ports its Bayesian
network inference, reconstruction, layout, visualization, and C++
extension interfaces directly to R with **no Python runtime
dependency**.

All exact model description lengths, negative log-likelihoods, and
incremental MCMC updates are numerically verified to machine precision
(`< 1e-11`) against a pinned `graph-tool` 3.9 reference build.

## Installation

`rgraphtool` requires R with a C++23-capable compiler, plus the `Rcpp`
and `BH` (`>= 1.81.0`) packages. Optional packages `igraph`, `ggplot2`,
and `svglite` enable `igraph` conversion,
[`ggplot_graph()`](https://diogro.github.io/rgraphtool/reference/plotting_integration.md),
and dependency-free SVG export.

``` r

# From a local checkout:
install.packages(".", repos = NULL, type = "source")

# Or from GitHub using remotes / pak:
# remotes::install_github("diogro/rgraphtool")
```

> **Note:** `rgraphtool` uses portable C++23 `std::thread` parallelism
> (no `libomp` required on macOS) and requires a C++23-capable
> toolchain.

## Quick start

### 1. Creating graphs & property maps

``` r

library(rgraphtool)

g <- graph_create(4L, directed = TRUE)
graph_add_edges(
  g,
  rbind(c(1, 2), c(2, 3), c(3, 1), c(2, 4)),
  properties = list(weight = c(1.2, 0.8, 2.5, 1.0), label = c("a", "b", "c", "d"))
)
graph_set_vertex_property(g, "name", c("v1", "v2", "v3", "v4"))
graph_set_graph_property(g, "title", "toy network")

g                              # summary header + compact ASCII adjacency for small graphs
graph_degree(g, "out")         # 1 2 1 0
graph_edge_property(g, "weight")
```

Vertex IDs are 1-based. Parallel edges and self-loops are supported.
Graphs are mutable external-pointer objects; persist graphs across
sessions using
[`write_gt()`](https://diogro.github.io/rgraphtool/reference/graph_io.md)
/
[`read_gt()`](https://diogro.github.io/rgraphtool/reference/graph_io.md)
(below).

### 2. Graph I/O: `graph-tool` (`.gt`) and `igraph` interoperability

`rgraphtool` implements the official `graph-tool` `.gt` binary format
(v1, little- and big-endian, all 15 property map types, with optional
`.gz`, `.bz2`, `.xz`, or `.zst` compression) natively in C++ with no
Python dependency, as well as lossless bidirectional conversion with
`igraph` (plus GraphML, GML, and edgelist files via
[`graph_load()`](https://diogro.github.io/rgraphtool/reference/graph_io.md)
/
[`graph_save()`](https://diogro.github.io/rgraphtool/reference/graph_io.md)):

``` r

# Native C++ graph-tool .gt binary read/write (interoperable with Python graph-tool)
write_gt(g, "network.gt.gz")
g2 <- read_gt("network.gt.gz")

# Bidirectional conversion with igraph (preserving graph, vertex, and edge attributes)
if (requireNamespace("igraph", quietly = TRUE)) {
  ig <- as_igraph(g)
  g_from_ig <- from_igraph(ig)
}

# Format-detecting load/save (.gt, .gt.gz, .graphml, .gml, .edgelist)
graph_save(g, "network.graphml")
g3 <- graph_load("network.graphml")
```

### 3. Hierarchical & bipartite SBM inference + posterior sampling

``` r

# Hierarchical SBM on C. elegans from Netzschleuder (https://networks.skewed.de)
celegans <- netzschleuder_network("celegansneural")

set.seed(42)
fit <- minimize_nested_blockmodel_dl(celegans$graph)
fit
plot_hierarchy(fit)

# Bipartite SBM directly from an ecological incidence matrix (154 fungi x 51 trees)
ft <- fungus_tree_network()
fit_bip <- minimize_bipartite_blockmodel_dl(ft$fungus_tree, dim_labels = c("fungus", "tree"))
plot_bipartite(fit_bip, type = "matrix")

# Posterior MCMC sampling and partition mode clustering
samples <- sample_posterior(fit, nsamples = 100L, burnin = 50L, thin = 5L)
modes <- mode_cluster_state(samples)
mode_summary(modes)
```

### 4. Network reconstruction from dynamics

``` r

# Reconstruct an interaction graph directly from equilibrium Ising spins
g_true <- graph_create(20L, directed = FALSE)
graph_add_edges(g_true, rbind(cbind(1:19, 2:20), c(20L, 1L)))

set.seed(1)
spins <- simulate_ising(g_true, weights = 0.6, n_samples = 400L)
st <- pseudo_ising_block_state(spins, nested = TRUE)
rec <- reconstruct_dynamics(st, niter = 20L)
rec
```

### 5. Layouts, animations, interactive views & maps

``` r

set.seed(7)
pos <- sfdp_layout(g_true)
plot_graph(g_true, pos = pos, edge_curve = 0.12, axes = TRUE,
           output = tempfile(fileext = ".pdf"))

# Standalone zero-dependency HTML5 Canvas animations & interactive BFS explorer
animate_sirs(g_true, pos = pos, n_steps = 60L, mode = "html", open = TRUE)
interactive_bfs(g_true, pos = pos, max_dist = 3L, mode = "html", open = TRUE)
```

## Capabilities & articles

The package website
([diogro.github.io/rgraphtool](https://diogro.github.io/rgraphtool/))
uses this README as its home page and hosts 9 grouped cookbook articles
covering all 22 bundled demos:

| Module | Key functions | Article |
|----|----|----|
| **Core graphs, I/O, flat & nested SBMs** | [`graph_create()`](https://diogro.github.io/rgraphtool/reference/graph.md), [`graph_add_edges()`](https://diogro.github.io/rgraphtool/reference/graph.md), [`graph_set_vertex_property()`](https://diogro.github.io/rgraphtool/reference/graph.md), [`graph_set_edge_property()`](https://diogro.github.io/rgraphtool/reference/graph.md), [`read_gt()`](https://diogro.github.io/rgraphtool/reference/graph_io.md), [`write_gt()`](https://diogro.github.io/rgraphtool/reference/graph_io.md), [`as_igraph()`](https://diogro.github.io/rgraphtool/reference/graph_io.md), [`from_igraph()`](https://diogro.github.io/rgraphtool/reference/graph_io.md), [`graph_load()`](https://diogro.github.io/rgraphtool/reference/graph_io.md), [`graph_save()`](https://diogro.github.io/rgraphtool/reference/graph_io.md), [`netzschleuder_network()`](https://diogro.github.io/rgraphtool/reference/netzschleuder_network.md), [`block_state()`](https://diogro.github.io/rgraphtool/reference/block_state.md), [`minimize_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/minimize_blockmodel_dl.md), [`compare_blockmodels()`](https://diogro.github.io/rgraphtool/reference/minimize_blockmodel_dl.md), [`nested_block_state()`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md), [`minimize_nested_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md), [`plot_hierarchy()`](https://diogro.github.io/rgraphtool/reference/plot_hierarchy.md) | [Block models & hierarchy](https://diogro.github.io/rgraphtool/articles/blockmodels-and-hierarchy.html) |
| **Posterior MCMC & partition modes** | [`mcmc_sweep()`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md), [`mcmc_equilibrate()`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md), [`mcmc_anneal()`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md), [`sample_posterior()`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md), [`align_partition_labels()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md), [`partition_overlap()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md), [`partition_mode_state()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md), [`mode_marginals()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md), [`mode_consensus()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md), [`mode_cluster_state()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md), [`plot_marginals()`](https://diogro.github.io/rgraphtool/reference/plot_marginals.md) | [Posterior sampling](https://diogro.github.io/rgraphtool/articles/posterior-sampling.html) |
| **Overlapping SBMs & model evidence** | [`overlap_block_state()`](https://diogro.github.io/rgraphtool/reference/overlap_block_state.md), [`minimize_overlap_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/overlap_block_state.md), [`nested_overlap_block_state()`](https://diogro.github.io/rgraphtool/reference/overlap_block_state.md), [`minimize_nested_overlap_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/overlap_block_state.md), [`overlap_mcmc_sweep()`](https://diogro.github.io/rgraphtool/reference/overlap_block_state.md), [`model_evidence()`](https://diogro.github.io/rgraphtool/reference/compare_model_classes.md), [`compare_model_classes()`](https://diogro.github.io/rgraphtool/reference/compare_model_classes.md) | [Model class selection](https://diogro.github.io/rgraphtool/articles/model-class-selection.html) |
| **Weighted & layered SBMs** | [`weighted_block_state()`](https://diogro.github.io/rgraphtool/reference/weighted_block_state.md), [`minimize_weighted_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/weighted_block_state.md), [`nested_weighted_block_state()`](https://diogro.github.io/rgraphtool/reference/weighted_block_state.md), [`transform_weights()`](https://diogro.github.io/rgraphtool/reference/graph.md), [`compare_weight_models()`](https://diogro.github.io/rgraphtool/reference/compare_weight_models.md), [`layered_block_state()`](https://diogro.github.io/rgraphtool/reference/layered_block_state.md), [`minimize_layered_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/layered_block_state.md), [`nested_layered_block_state()`](https://diogro.github.io/rgraphtool/reference/layered_block_state.md), [`compare_layered_models()`](https://diogro.github.io/rgraphtool/reference/compare_layered_models.md) | [Weighted & layered SBMs](https://diogro.github.io/rgraphtool/articles/weighted-and-layered-sbms.html) |
| **Assortative & ordered (ranked) SBMs** | [`pp_block_state()`](https://diogro.github.io/rgraphtool/reference/pp_block_state.md), [`minimize_pp_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/pp_block_state.md), [`pp_assortativity()`](https://diogro.github.io/rgraphtool/reference/pp_block_state.md), [`compare_assortative_models()`](https://diogro.github.io/rgraphtool/reference/pp_block_state.md), [`ranked_block_state()`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md), [`minimize_ranked_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md), [`nested_ranked_block_state()`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md), [`ranked_vertex_order()`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md), [`ranked_edge_directions()`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md), [`compare_ordered_models()`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md), [`plot_ranked_hierarchy()`](https://diogro.github.io/rgraphtool/reference/plot_hierarchy.md) | [Assortative & ordered SBMs](https://diogro.github.io/rgraphtool/articles/assortative-and-ordered-sbms.html) |
| **Bipartite & ecological SBMs** | [`bipartite_graph()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md), [`is_bipartite()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md), [`as_incidence_matrix()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md), [`fungus_tree_network()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md), [`bipartite_block_state()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md), [`minimize_bipartite_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md), [`nested_bipartite_block_state()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md), [`minimize_nested_bipartite_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md), [`bipartite_membership()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md), [`bipartite_block_matrix()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md), [`compare_bipartite_models()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md), [`plot_bipartite()`](https://diogro.github.io/rgraphtool/reference/plot_bipartite.md) | [Bipartite & ecological networks](https://diogro.github.io/rgraphtool/articles/bipartite-sbms.html) |
| **Uncertain networks, triadic closure & dynamics** | [`measured_block_state()`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md), [`mixed_measured_block_state()`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md), [`uncertain_block_state()`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md), [`latent_multigraph()`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md), [`reconstruct_network()`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md), [`predict_edges()`](https://diogro.github.io/rgraphtool/reference/predict_edges.md), [`latent_closure_block_state()`](https://diogro.github.io/rgraphtool/reference/latent_closure_block_state.md), [`infer_triadic_closure()`](https://diogro.github.io/rgraphtool/reference/latent_closure_block_state.md), all 10 dynamical models (`pseudo_ising`, `ising_glauber`, `pseudo_cising`, `cising_glauber`, `epidemics`, `pseudo_normal`, `normal_glauber`, `linear_normal`, `pseudo_potts`, `potts_glauber`), [`reconstruct_dynamics()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md), [`reconstruct_l1_path()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md), [`simulate_ising()`](https://diogro.github.io/rgraphtool/reference/simulate_dynamics.md), [`simulate_epidemics()`](https://diogro.github.io/rgraphtool/reference/simulate_dynamics.md), [`simulate_normal_dynamics()`](https://diogro.github.io/rgraphtool/reference/simulate_dynamics.md), [`simulate_potts()`](https://diogro.github.io/rgraphtool/reference/simulate_dynamics.md) | [Network reconstruction](https://diogro.github.io/rgraphtool/articles/network-reconstruction.html) |
| **Layouts, animations, plotting & maps** | [`sfdp_layout()`](https://diogro.github.io/rgraphtool/reference/sfdp_layout.md), [`graph_bfs()`](https://diogro.github.io/rgraphtool/reference/sfdp_layout.md), [`plot_graph()`](https://diogro.github.io/rgraphtool/reference/plot_graph.md), [`animate_sirs()`](https://diogro.github.io/rgraphtool/reference/animate_sirs.md), [`animate_dynamic_layout()`](https://diogro.github.io/rgraphtool/reference/animate_sirs.md), [`interactive_bfs()`](https://diogro.github.io/rgraphtool/reference/interactive_bfs.md), [`prop_to_size()`](https://diogro.github.io/rgraphtool/reference/plotting_integration.md), [`map_colors()`](https://diogro.github.io/rgraphtool/reference/plotting_integration.md), [`plot_colorbar()`](https://diogro.github.io/rgraphtool/reference/plotting_integration.md), [`save_graph_plot()`](https://diogro.github.io/rgraphtool/reference/plotting_integration.md), [`graph_grob()`](https://diogro.github.io/rgraphtool/reference/plotting_integration.md), [`ggplot_graph()`](https://diogro.github.io/rgraphtool/reference/plotting_integration.md), [`project_coords()`](https://diogro.github.io/rgraphtool/reference/plot_graph_map.md), [`great_circle_arc()`](https://diogro.github.io/rgraphtool/reference/plot_graph_map.md), [`world_coastlines()`](https://diogro.github.io/rgraphtool/reference/plot_graph_map.md), [`plot_graph_map()`](https://diogro.github.io/rgraphtool/reference/plot_graph_map.md) | [Visualization & animation](https://diogro.github.io/rgraphtool/articles/visualization-and-animation.html) |
| **Large-network scaling & C++ extensions** | [`rgraphtool_max_threads()`](https://diogro.github.io/rgraphtool/reference/rgraphtool_threads.md), [`rgraphtool_num_threads()`](https://diogro.github.io/rgraphtool/reference/rgraphtool_threads.md), [`rgraphtool_set_num_threads()`](https://diogro.github.io/rgraphtool/reference/rgraphtool_threads.md), public `<rgraphtool.h>` API (`LinkingTo: rgraphtool, Rcpp, BH`), [`kcore_decomposition()`](https://diogro.github.io/rgraphtool/reference/cpp_extensions.md), [`vertex_strength()`](https://diogro.github.io/rgraphtool/reference/cpp_extensions.md), [`rgraphtool_include_dir()`](https://diogro.github.io/rgraphtool/reference/cpp_extensions.md), [`rgraphtool_extension_example()`](https://diogro.github.io/rgraphtool/reference/cpp_extensions.md), [`compile_cpp_extension()`](https://diogro.github.io/rgraphtool/reference/cpp_extensions.md) | [Scaling & C++ extensions](https://diogro.github.io/rgraphtool/articles/scaling-and-cpp-extensions.html) |

All 22 underlying cookbook scripts remain runnable interactively via
`demo("<name>", package = "rgraphtool")` (see
`demo(package = "rgraphtool")`).

## Multi-threading for large networks

`rgraphtool` includes a portable header-only C++23 thread pool
(`src/parallel.h`) that works across macOS, Linux, and Windows without
requiring OpenMP:

``` r

rgraphtool_max_threads()
rgraphtool_set_num_threads(4L)
# Or via global option:
options(rgraphtool.nthreads = 4L)
```

## Writing C++ extensions

User C++ extensions and downstream R packages
(`LinkingTo: rgraphtool, Rcpp, BH`) can `#include <rgraphtool.h>` for
zero-copy access to `rgraphtool::Graph`, range-based graph iteration
(`vertices_range`, `edges_range`, `out_edges_range`, `in_edges_range`,
`all_edges_range`, `out_neighbors_range`, `in_neighbors_range`,
`all_neighbors_range`), and vertex/edge property maps with
bounds-checked defaults and upfront-validated zero-overhead unchecked
access (`.get_unchecked()`):

``` r

# Compile and run a bundled standalone C++ extension via Rcpp::sourceCpp()
compile_cpp_extension(rgraphtool_extension_example("kcore"))
ext_kcore_decomposition(g, checked = FALSE)
```

See
[`?cpp_extensions`](https://diogro.github.io/rgraphtool/reference/cpp_extensions.md)
and `demo("cpp_extensions", package = "rgraphtool")` for complete
examples and benchmarks.

## Validation against graph-tool 3.9

`tools/reference/` contains a pinned `pixi` environment with the
conda-forge build of `graph-tool` 3.9
(`b24bc428e8426d0acb36b2d57bd61633335a3a26`) along with the scripts used
to generate numerical reference fixtures for the test suite (1,272
tests). Python is **never** required to install, use, or test
`rgraphtool`.

## Attribution and license

`rgraphtool` is derived from [graph-tool](https://graph-tool.skewed.de/)
release 3.9 by Tiago P. Peixoto. Vendored adjacency headers reside in
`inst/include/rgraphtool/vendor/` with their original copyright notices
intact (unused generated `config.h` includes were removed from two
supporting headers). Both `graph-tool` and `rgraphtool` are distributed
under the GNU Lesser General Public License, version 3 or later
(`LGPL-3.0-or-later`); see `COPYING` and `COPYING.LESSER`.
