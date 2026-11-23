# rgraphtool

`rgraphtool` is an experimental, incremental **native R/C++ port** of the Python
[graph-tool](https://graph-tool.skewed.de/) library (originally developed by
Tiago P. Peixoto). It incorporates `graph-tool`'s C++ adjacency structure and
ports its Bayesian stochastic block model (SBM) inference engines into native
C++ for R, with **no Python runtime dependency**.

Because `rgraphtool` translates `graph-tool`'s C++ data structures and
algorithms to R's native ecosystem (rather than wrapping Python via
`reticulate` or maintaining a Python fork), it is best described as a native
R/C++ **port** of Python `graph-tool`. It is **not** yet a complete port of
all `graph-tool` submodules (see `ROADMAP.md` for completed and planned
milestones, including large-network scaling and uncertain network
reconstruction).

The active adjacency headers come from `graph-tool` release 3.9, commit
`b24bc428e8426d0acb36b2d57bd61633335a3a26`. Their original copyright
notices are retained. The unused generated `config.h` includes were removed
from two supporting headers. `graph-tool` and this project are distributed
under LGPL-3.0-or-later; see `COPYING` and `COPYING.LESSER`.

## Installation

Install R with a configured C++23 toolchain, and the R packages `Rcpp` and
`BH` (version 1.81.0 or later), then run:

```r
install.packages(".", repos = NULL, type = "source")
```

Not all R platforms currently support C++23.

## Example

```r
library(rgraphtool)
g <- graph_create(3)
graph_add_edges(g, rbind(c(1, 2), c(2, 3), c(2, 3)))
graph_num_edges(g)             # 3
graph_degree(g, "out")        # 1 2 0
graph_edges(g)
```

Vertex IDs are 1-based. Parallel edges and self-loops are supported.
Undirected graphs treat stored endpoints as unordered for degree purposes;
`graph_edges()` returns them in their insertion orientation. Total degree counts
a self-loop twice. Graphs are mutable: aliases refer to the same graph.
External pointers are process-local and cannot be saved and restored with R
serialization. Keep the edges separately if persistence is needed.

## Implemented inference models

All description lengths are verified to machine precision (`< 1e-11`) against
pinned `graph-tool` 3.9 reference outputs:

- **Flat and nested (hierarchical) SBMs**: `block_state()`,
  `minimize_blockmodel_dl()`, `nested_block_state()`,
  `minimize_nested_blockmodel_dl()`, `compare_blockmodels()`,
  `plot_hierarchy()`. See `demo("football")`, `demo("celegans")`,
  `demo("model_selection")`.
- **MCMC posterior sampling & partition modes**: `mcmc_sweep()`,
  `mcmc_equilibrate()`, `mcmc_anneal()`, `sample_posterior()`,
  `partition_mode_state()`, `mode_cluster_state()`, `plot_marginals()`. See
  `demo("posterior_sampling")`, `demo("hierarchical_sampling")`.
- **Overlapping SBMs & Bayesian model class selection**:
  `overlap_block_state()`, `minimize_overlap_blockmodel_dl()`,
  `nested_overlap_block_state()`, `minimize_nested_overlap_blockmodel_dl()`,
  `model_evidence()`, `compare_model_classes()`. See
  `demo("model_class_selection")`.
- **Weighted / edge-covariate SBMs**: edge property maps
  (`graph_set_edge_property()`, `graph_edge_property()`),
  `weighted_block_state()`, `minimize_weighted_blockmodel_dl()`,
  `nested_weighted_block_state()`, `minimize_nested_weighted_blockmodel_dl()`,
  `transform_weights()`, `compare_weight_models()`. See `demo("weighted_sbm")`.
- **Layered SBMs**: `layered_block_state()`, `minimize_layered_blockmodel_dl()`,
  `nested_layered_block_state()`, `minimize_nested_layered_blockmodel_dl()`,
  `compare_layered_models()`. See `demo("layered_sbm")`.
- **Assortative (planted-partition) & ordered (ranked) SBMs**:
  `pp_block_state()`, `minimize_pp_blockmodel_dl()`,
  `compare_assortative_models()`, `ranked_block_state()`,
  `minimize_ranked_blockmodel_dl()`, `nested_ranked_block_state()`,
  `minimize_nested_ranked_blockmodel_dl()`, `compare_ordered_models()`,
  `plot_ranked_hierarchy()`. See `demo("assortative_sbm")`,
  `demo("ordered_sbm")`.

`netzschleuder_network()` downloads any network from the
[Netzschleuder](https://networks.skewed.de) repository used by `graph-tool`'s
`collection.ns`.

```r
library(rgraphtool)
g <- graph_create(3, directed = FALSE)
graph_add_edges(g, rbind(c(1, 2), c(2, 3)))
state <- block_state(g, c(1, 1, 2))
block_edge_counts(state)
block_log_likelihood(state)
block_description_length(state)

football <- football_network()
set.seed(1)
fit <- minimize_blockmodel_dl(football$graph)
block_sizes(fit)
block_description_length(fit)

celegans <- netzschleuder_network("celegansneural")
nested <- minimize_nested_blockmodel_dl(celegans$graph)
nested
plot_hierarchy(nested)
```

## Validation against graph-tool

`tools/reference/` holds a pixi environment pinned to the conda-forge build of
`graph-tool` 3.9, plus the scripts that produced the reference values in the
tests. Python is only needed to regenerate them; see
`tools/reference/README.md`.
