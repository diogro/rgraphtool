# Inspect a fixed block-model partition

Build and inspect a fixed block-model partition and score it with the
microcanonical stochastic block model description length of graph-tool
3.9. See
[`minimize_blockmodel_dl`](https://diogro.github.io/rgraphtool/reference/minimize_blockmodel_dl.md)
to fit a partition.

## Usage

``` r
block_state(graph, blocks = NULL, degree_corrected = TRUE)
block_membership(state)
block_sizes(state)
block_edge_counts(state)
block_log_likelihood(state)
block_description_length(state, degree_prior = NULL)
```

## Arguments

- graph:

  A graph created by
  [`graph_create()`](https://diogro.github.io/rgraphtool/reference/graph.md).

- blocks:

  One positive integer group label per vertex; defaults to a single
  group.

- degree_corrected:

  Whether to use the degree-corrected likelihood.

- state:

  A state created by `block_state()`.

- degree_prior:

  Degree-sequence prior, either `"distributed"` or `"uniform"`. For
  `block_description_length()`, `NULL` uses the prior associated with
  the state.

## Details

Matrix rows and columns correspond to the sorted nonempty labels.
Directed counts are source by target. For undirected graphs the matrix
is symmetric, and diagonal values count each within-group edge twice,
including self-loops. Parallel edges and self-loops are included. The
likelihood is the log probability of the observed multigraph conditional
on the block counts, partition, and, for degree correction, degree
sequence. It does not include priors or description length.
`block_description_length()` adds the partition, edge-count, and, for
degree-corrected states, degree-sequence priors, matching graph-tool's
`BlockState.entropy()` with its default arguments (`dense = FALSE`,
`multigraph = TRUE`, `exact = TRUE`). The partition prior includes the
uniform prior on the number of groups, \\\log N\\.

The distributed degree prior needs the number of integer partitions
\\q(m, n)\\. As in graph-tool, exact values are tabulated up to
\\\min(\max(2E, 100), 10000)\\ and an asymptotic approximation is used
above that. The exact table needs about \\4 n^2\\ bytes, so this package
caps it at `getOption("rgraphtool.q_cache_max", 5000)` (about 100 MB);
raise the option to 10000 for exact agreement with graph-tool whenever a
group's degree sum lies between 5000 and 10000. The table is kept for
the session.

States are snapshots: later graph changes do not change their
statistics. If the graph has grown since the state was created, the
description length is computed from the stored snapshot in R (exact
partition counts, limited to 10 million table entries) instead of the
native code.

## Value

`block_state()` returns a block state. `block_membership()` returns
vertex labels, `block_sizes()` returns counts by label,
`block_edge_counts()` returns the affinity matrix, and
`block_log_likelihood()` returns a numeric scalar in nats.
`block_description_length()` returns the description length in nats.

## See also

[`minimize_blockmodel_dl`](https://diogro.github.io/rgraphtool/reference/minimize_blockmodel_dl.md)

## Examples

``` r
g <- graph_create(3, directed = FALSE)
graph_add_edges(g, rbind(c(1, 2), c(2, 3)))
state <- block_state(g, c(1, 1, 2))
block_edge_counts(state)
#>   1 2
#> 1 2 1
#> 2 1 0
block_log_likelihood(state)
#> [1] -0.4054651
block_description_length(state)
#> [1] 6.473891
```
