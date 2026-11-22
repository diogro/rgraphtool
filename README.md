# graph-tool-R

An experimental, incremental native R port of graph-tool. It uses graph-tool's C++
adjacency structure and implements flat stochastic block model inference; it is
**not** a complete port and does not yet implement property maps, filters,
hierarchical inference, file formats, or visualization. It does not depend on Python.

The active adjacency headers come from graph-tool release 3.9, commit
`b24bc428e8426d0acb36b2d57bd61633335a3a26`. Their original copyright
notices are retained. The unused generated `config.h` includes were removed
from two supporting headers.
graph-tool and this project are distributed under LGPL-3.0-or-later; see
`COPYING` and `COPYING.LESSER`.

## Installation

Install R with a configured C++23 toolchain, and the R packages `Rcpp` and
`BH` (version 1.81.0 or later), then run:

```r
install.packages(".", repos = NULL, type = "source")
```

Not all R platforms currently support C++23.

## Example

```r
library(graphtool)
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

## Stochastic block models

`block_state()` records a supplied partition (or a single group by default),
with block sizes, edge counts, conditional microcanonical SBM log likelihood,
and description length. Description lengths are computed natively and match
graph-tool 3.9's `BlockState.entropy()` to a relative tolerance of 1e-10, for
directed and undirected multigraphs, with and without degree correction.

`minimize_blockmodel_dl()` fits a flat SBM with a native port of graph-tool's
multilevel agglomerative search: vertex sweeps from singleton groups, a
bisection search over the number of groups with merges and zero-temperature
sweeps, then refinement. On the cookbook football network it gives the same
distribution of description lengths as graph-tool 3.9 (best 1738.87 nats
non-degree-corrected, 1787.36 degree-corrected) in a few hundredths of a
second. Results are local optima: fit several seeds and keep the smallest
description length. `refine_blockmodel_dl()` runs only the refinement stage on
a supplied partition, and `compare_blockmodels()` fits both models and reports
the one with the smaller description length. See `?minimize_blockmodel_dl` for
the algorithm and its differences from graph-tool.

Nested (hierarchical) models, posterior sampling, and other cookbook datasets
are not implemented yet. The football network loader and base-R fit/plot
example are available via `football_network()` and `demo("football")`.

```r
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
```

## Validation against graph-tool

`tools/reference/` holds a pixi environment pinned to the conda-forge build of
graph-tool 3.9, plus the scripts that produced the reference values in the
tests. Python is only needed to regenerate them; see
`tools/reference/README.md`.
