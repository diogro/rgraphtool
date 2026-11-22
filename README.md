# graph-tool-R

An experimental, incremental native R port of graph-tool. This first milestone
uses graph-tool's C++ adjacency structure; it is **not** a complete port and
does not yet implement property maps, filters, graph algorithms, file formats,
or visualization. It does not depend on Python.

The vendored adjacency implementation comes from the graph-tool checkout at
commit `85f11ae869fe318532e60554963b16eca6367eab` (2015-06-18), not the
current graph-tool release. Its original copyright notices are retained.
The headers were renamed for R package conventions. graph-tool and this project
are distributed under GPL-3.0-or-later; see
`COPYING`.

## Installation

Install R, a C++17 toolchain, and the R packages `Rcpp` and `BH`, then run:

```r
install.packages(".", repos = NULL, type = "source")
```

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
