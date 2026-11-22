# graph-tool-R

An experimental, incremental native R port of graph-tool. This first milestone
uses graph-tool's C++ adjacency structure; it is **not** a complete port and
does not yet implement property maps, filters, hierarchical inference, file formats,
or visualization. It does not depend on Python.

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

## Block-model foundations

`block_state()` records a supplied partition (or a single group by default),
with block sizes, edge counts, conditional microcanonical SBM log likelihood,
and description length. `minimize_blockmodel_dl()` performs an experimental
greedy flat-SBM fit. It is not graph-tool's multilevel optimizer, may find
worse local optima, and is intended for small graphs. See `?block_state`.
`compare_blockmodels()` fits both degree-corrected and non-degree-corrected
flat models and reports which has the smaller description length.
`refine_blockmodel_dl()` greedily refines a supplied flat partition without
increasing its description length; it is not graph-tool's merge-split MCMC.
Hierarchical inference, posterior sampling, and other cookbook datasets are
not implemented yet. The football network loader and base-R fit/plot example
are available via `football_network()` and `demo("football")`.

```r
g <- graph_create(3, directed = FALSE)
graph_add_edges(g, rbind(c(1, 2), c(2, 3)))
state <- block_state(g, c(1, 1, 2))
block_edge_counts(state)
block_log_likelihood(state)
block_description_length(state)
set.seed(1)
fit <- minimize_blockmodel_dl(g)
block_membership(fit)
```
