# Draw a network with pie-chart vertex marginal group probabilities

Base-R visualization of posterior vertex marginal group memberships,
drawing a pie chart at each vertex proportional to its marginal
probabilities across aligned groups (the R equivalent of graph-tool's
`vertex_shape="pie"` marginal plots).

## Usage

``` r
plot_marginals(graph, mode, positions = NULL, radius = 0.04,
  palette = NULL, edge_col = grDevices::adjustcolor("grey40", 0.35),
  border = "grey20", labels = NULL, label_cex = 0.65, main = NULL, ...)
```

## Arguments

- graph:

  A graph created by
  [`graph_create`](https://diogro.github.io/rgraphtool/reference/graph.md).

- mode:

  A `rgraphtool_partition_mode` from
  [`partition_mode_state`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md)
  or
  [`mode_cluster_state`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md),
  or an \\N \times B\\ numeric matrix of marginal counts or
  probabilities.

- positions:

  Optional \\N \times 2\\ matrix of vertex coordinates. When `NULL`,
  computes a deterministic spring layout initialized from the consensus
  partition.

- radius:

  Pie-chart radius in normalized plot coordinates.

- palette:

  Optional vector of group colours.

- edge_col:

  Edge colour.

- border:

  Pie slice border colour.

- labels:

  Optional character vector of vertex labels.

- label_cex:

  Character expansion for `labels`.

- main:

  Optional plot title.

- ...:

  Additional arguments passed to
  [`graphics::plot.window()`](https://rdrr.io/r/graphics/plot.window.html).

## Value

Invisibly, a list with normalized vertex `positions`, the marginal
`probabilities` matrix, and the group colour `palette`.

## See also

[`partition_mode_state`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md),
[`plot_hierarchy`](https://diogro.github.io/rgraphtool/reference/plot_hierarchy.md),
`demo("posterior_sampling")`

## Examples

``` r
g <- graph_create(6, directed = FALSE)
graph_add_edges(g, rbind(c(1, 2), c(2, 3), c(1, 3), c(4, 5), c(5, 6), c(4, 6), c(3, 4)))
pm <- partition_mode_state(list(c(1, 1, 1, 2, 2, 2), c(1, 1, 2, 2, 2, 2)))
plot_marginals(g, pm, main = "Vertex marginals")
```
