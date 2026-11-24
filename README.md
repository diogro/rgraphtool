# rgraphtool

`rgraphtool` is a standalone **native R/C++ port** of the Python
[graph-tool](https://graph-tool.skewed.de/) network analysis library
(originally developed by Tiago P. Peixoto). It embeds `graph-tool`'s C++
adjacency data structure (`boost::adj_list`) and ports its Bayesian network
inference, reconstruction, layout, visualization, and C++ extension interfaces
directly to R with **no Python runtime dependency**.

All exact model description lengths, negative log-likelihoods, and incremental
MCMC updates are numerically verified to machine precision (`< 1e-11`) against
a pinned `graph-tool` 3.9 reference build.

## Installation

`rgraphtool` requires R with a C++23-capable compiler, plus the `Rcpp` and `BH`
(`>= 1.81.0`) packages. Optional packages `ggplot2` and `svglite` enable
`ggplot_graph()` and dependency-free SVG export.

```r
# From a local checkout:
install.packages(".", repos = NULL, type = "source")

# Or from GitHub using remotes / pak:
# remotes::install_github("diogro/rgraphtool")
```

> **Note:** `rgraphtool` uses portable C++23 `std::thread` parallelism (no
> `libomp` required on macOS) and requires a C++23-capable toolchain.

## Quick start

### 1. Creating graphs & edge properties

```r
library(rgraphtool)

g <- graph_create(4L, directed = TRUE)
graph_add_edges(
  g,
  rbind(c(1, 2), c(2, 3), c(3, 1), c(2, 4)),
  properties = list(weight = c(1.2, 0.8, 2.5, 1.0))
)
graph_num_vertices(g)          # 4
graph_num_edges(g)             # 4
graph_degree(g, "out")         # 1 2 1 0
graph_edge_property(g, "weight")
```

Vertex IDs are 1-based. Parallel edges and self-loops are supported. Graphs are
mutable external-pointer objects (process-local; store edge tables separately
if saving to disk with `saveRDS()`).

### 2. Hierarchical SBM inference & posterior sampling

```r
# Load any network from Netzschleuder (https://networks.skewed.de)
celegans <- netzschleuder_network("celegansneural")

set.seed(42)
fit <- minimize_nested_blockmodel_dl(celegans$graph)
fit
plot_hierarchy(fit)

# Posterior MCMC sampling and partition mode clustering
samples <- sample_posterior(fit, nsamples = 100L, burnin = 50L, thin = 5L)
modes <- mode_cluster_state(samples)
mode_summary(modes)
```

### 3. Network reconstruction from dynamics

```r
# Reconstruct an interaction graph directly from equilibrium Ising spins
g_true <- graph_create(20L, directed = FALSE)
graph_add_edges(g_true, rbind(cbind(1:19, 2:20), c(20L, 1L)))

set.seed(1)
spins <- simulate_ising(g_true, weights = 0.6, n_samples = 400L)
st <- pseudo_ising_block_state(spins, nested = TRUE)
rec <- reconstruct_dynamics(st, niter = 20L)
rec
```

### 4. Layouts, animations, interactive views & maps

```r
set.seed(7)
pos <- sfdp_layout(g_true)
plot_graph(g_true, pos = pos, edge_curve = 0.12, axes = TRUE,
           output = tempfile(fileext = ".pdf"))

# Standalone zero-dependency HTML5 Canvas animations & interactive BFS explorer
animate_sirs(g_true, pos = pos, n_steps = 60L, mode = "html", open = TRUE)
interactive_bfs(g_true, pos = pos, max_dist = 3L, mode = "html", open = TRUE)
```

## Capabilities & articles

The package website ([diogro.github.io/rgraphtool](https://diogro.github.io/rgraphtool/))
uses this README as its home page and hosts 8 grouped cookbook articles
covering all 21 bundled demos:

| Module | Key functions | Article |
| --- | --- | --- |
| **Core graphs, flat & nested SBMs** | `graph_create()`, `graph_add_edges()`, `graph_set_edge_property()`, `netzschleuder_network()`, `block_state()`, `minimize_blockmodel_dl()`, `compare_blockmodels()`, `nested_block_state()`, `minimize_nested_blockmodel_dl()`, `plot_hierarchy()` | [Block models & hierarchy](https://diogro.github.io/rgraphtool/articles/blockmodels-and-hierarchy.html) |
| **Posterior MCMC & partition modes** | `mcmc_sweep()`, `mcmc_equilibrate()`, `mcmc_anneal()`, `sample_posterior()`, `align_partition_labels()`, `partition_overlap()`, `partition_mode_state()`, `mode_marginals()`, `mode_consensus()`, `mode_cluster_state()`, `plot_marginals()` | [Posterior sampling](https://diogro.github.io/rgraphtool/articles/posterior-sampling.html) |
| **Overlapping SBMs & model evidence** | `overlap_block_state()`, `minimize_overlap_blockmodel_dl()`, `nested_overlap_block_state()`, `minimize_nested_overlap_blockmodel_dl()`, `overlap_mcmc_sweep()`, `model_evidence()`, `compare_model_classes()` | [Model class selection](https://diogro.github.io/rgraphtool/articles/model-class-selection.html) |
| **Weighted & layered SBMs** | `weighted_block_state()`, `minimize_weighted_blockmodel_dl()`, `nested_weighted_block_state()`, `transform_weights()`, `compare_weight_models()`, `layered_block_state()`, `minimize_layered_blockmodel_dl()`, `nested_layered_block_state()`, `compare_layered_models()` | [Weighted & layered SBMs](https://diogro.github.io/rgraphtool/articles/weighted-and-layered-sbms.html) |
| **Assortative & ordered (ranked) SBMs** | `pp_block_state()`, `minimize_pp_blockmodel_dl()`, `pp_assortativity()`, `compare_assortative_models()`, `ranked_block_state()`, `minimize_ranked_blockmodel_dl()`, `nested_ranked_block_state()`, `ranked_vertex_order()`, `ranked_edge_directions()`, `compare_ordered_models()`, `plot_ranked_hierarchy()` | [Assortative & ordered SBMs](https://diogro.github.io/rgraphtool/articles/assortative-and-ordered-sbms.html) |
| **Uncertain networks, triadic closure & dynamics** | `measured_block_state()`, `mixed_measured_block_state()`, `uncertain_block_state()`, `latent_multigraph()`, `reconstruct_network()`, `predict_edges()`, `latent_closure_block_state()`, `infer_triadic_closure()`, all 10 dynamical models (`pseudo_ising`, `ising_glauber`, `pseudo_cising`, `cising_glauber`, `epidemics`, `pseudo_normal`, `normal_glauber`, `linear_normal`, `pseudo_potts`, `potts_glauber`), `reconstruct_dynamics()`, `reconstruct_l1_path()`, `simulate_ising()`, `simulate_epidemics()`, `simulate_normal_dynamics()`, `simulate_potts()` | [Network reconstruction](https://diogro.github.io/rgraphtool/articles/network-reconstruction.html) |
| **Layouts, animations, plotting & maps** | `sfdp_layout()`, `graph_bfs()`, `plot_graph()`, `animate_sirs()`, `animate_dynamic_layout()`, `interactive_bfs()`, `prop_to_size()`, `map_colors()`, `plot_colorbar()`, `save_graph_plot()`, `graph_grob()`, `ggplot_graph()`, `project_coords()`, `great_circle_arc()`, `world_coastlines()`, `plot_graph_map()` | [Visualization & animation](https://diogro.github.io/rgraphtool/articles/visualization-and-animation.html) |
| **Large-network scaling & C++ extensions** | `rgraphtool_max_threads()`, `rgraphtool_num_threads()`, `rgraphtool_set_num_threads()`, public `<rgraphtool.h>` API (`LinkingTo: rgraphtool, Rcpp, BH`), `kcore_decomposition()`, `vertex_strength()`, `rgraphtool_include_dir()`, `rgraphtool_extension_example()`, `compile_cpp_extension()` | [Scaling & C++ extensions](https://diogro.github.io/rgraphtool/articles/scaling-and-cpp-extensions.html) |

All 21 underlying cookbook scripts remain runnable interactively via
`demo("<name>", package = "rgraphtool")` (see `demo(package = "rgraphtool")`).

## Multi-threading for large networks

`rgraphtool` includes a portable header-only C++23 thread pool (`src/parallel.h`)
that works across macOS, Linux, and Windows without requiring OpenMP:

```r
rgraphtool_max_threads()
rgraphtool_set_num_threads(4L)
# Or via global option:
options(rgraphtool.nthreads = 4L)
```

## Writing C++ extensions

User C++ extensions and downstream R packages (`LinkingTo: rgraphtool, Rcpp, BH`)
can `#include <rgraphtool.h>` for zero-copy access to `rgraphtool::Graph`,
range-based graph iteration (`vertices_range`, `edges_range`, `out_edges_range`,
`in_edges_range`, `all_edges_range`, `out_neighbors_range`, `in_neighbors_range`,
`all_neighbors_range`), and vertex/edge property maps with bounds-checked
defaults and upfront-validated zero-overhead unchecked access (`.get_unchecked()`):

```r
# Compile and run a bundled standalone C++ extension via Rcpp::sourceCpp()
compile_cpp_extension(rgraphtool_extension_example("kcore"))
ext_kcore_decomposition(g, checked = FALSE)
```

See `?cpp_extensions` and `demo("cpp_extensions", package = "rgraphtool")` for
complete examples and benchmarks.

## Validation against graph-tool 3.9

`tools/reference/` contains a pinned `pixi` environment with the conda-forge
build of `graph-tool` 3.9 (`b24bc428e8426d0acb36b2d57bd61633335a3a26`) along with
the scripts used to generate numerical reference fixtures for the test suite
(1,272 tests). Python is **never** required to install, use, or test `rgraphtool`.

## Attribution and license

`rgraphtool` is derived from [graph-tool](https://graph-tool.skewed.de/) release
3.9 by Tiago P. Peixoto. Vendored adjacency headers reside in
`inst/include/rgraphtool/vendor/` with their original copyright notices intact
(unused generated `config.h` includes were removed from two supporting headers).
Both `graph-tool` and `rgraphtool` are distributed under the GNU Lesser General
Public License, version 3 or later (`LGPL-3.0-or-later`); see `COPYING` and
`COPYING.LESSER`.
