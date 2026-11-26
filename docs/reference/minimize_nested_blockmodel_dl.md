# Hierarchical (nested) stochastic block models

The nested SBM places the groups of each level into groups at the level
above, up to a single group. Each upper level is a non-degree-corrected
SBM of the multigraph of group-pair edge counts below it. These
functions are the R counterparts of graph-tool's `NestedBlockState` and
`minimize_nested_blockmodel_dl()`.

## Usage

``` r
minimize_nested_blockmodel_dl(graph, degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"), B_min = 1, B_max = Inf,
  refine = TRUE, epsilon = 1e-3, simple_init = FALSE,
  simple_init_reverse = FALSE, base_factor = 10, B_min_base = 1000L,
  top_factor = 2, strict = FALSE, init_levels = NULL, nthreads = NULL,
  control = list())
nested_block_state(graph, levels, degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"))
nested_description_length(state, degree_prior = NULL)
nested_level_entropy(state, degree_prior = NULL)
nested_levels(state)
nested_project(state, level = 0L)
nested_level_state(state, level = 0L)
nested_summary(state)
```

## Arguments

- graph:

  A graph created by
  [`graph_create()`](https://diogro.github.io/rgraphtool/reference/graph.md),
  with at least one vertex.

- degree_corrected:

  Whether the lowest level uses the degree-corrected model.

- degree_prior:

  Degree-sequence prior of the lowest level. `NULL` uses the prior
  stored in `state`.

- B_min, B_max:

  Bounds on the number of groups at level 0.

- refine:

  Whether to run the refinement loop after the level-by-level search.

- epsilon:

  Relative change of the description length that stops refinement.

- simple_init:

  If `TRUE`, use graph-tool's fast hierarchy-shape initialization for
  large networks (strictly agglomerative level-by-level coarsening
  targeting `base_factor`, `B_min_base`, and `top_factor`) instead of
  full nested bisection search.

- simple_init_reverse:

  If `TRUE` with `simple_init = TRUE`, sweep levels top-down rather than
  bottom-up.

- base_factor:

  When `simple_init = TRUE`, target level-0 group count
  `max(N / base_factor, B_min_base)`.

- B_min_base:

  Minimum target level-0 group count during `simple_init = TRUE`
  initialization.

- top_factor:

  When `simple_init = TRUE`, shrink each upper level by
  `N_l / top_factor` (must be greater than 1).

- strict:

  When `simple_init = TRUE`, treat the target group counts as exact
  targets rather than lower bounds during initialization.

- init_levels:

  Optional warm-start hierarchy (either a list of 1-based label vectors
  or an existing `rgraphtool_nested_block_state`), skipping
  initialization and proceeding directly to refinement.

- nthreads:

  Optional positive integer number of parallel worker threads (or
  `"auto"`). `NULL` uses
  [`rgraphtool_num_threads`](https://diogro.github.io/rgraphtool/reference/rgraphtool_threads.md).

- control:

  Named list of search settings, as in
  [`minimize_blockmodel_dl`](https://diogro.github.io/rgraphtool/reference/minimize_blockmodel_dl.md).
  Here `c` defaults to 1 and is scaled by \\(2^l\\) at level \\(l\\), as
  in graph-tool.

- levels:

  A list of label vectors. `levels[[1]]` labels each vertex;
  `levels[[l]]` labels each group of `levels[[l - 1]]`, so it has one
  entry per distinct label of the previous level, in order of first
  appearance.

- state:

  A nested block state.

- level:

  Hierarchy level, counted from 0 (the vertex partition) as in
  graph-tool.

## Details

The description length is that of graph-tool 3.9's
`NestedBlockState.entropy()` with default arguments. Level 0 is the flat
SBM of
[`block_state`](https://diogro.github.io/rgraphtool/reference/block_state.md)
without the edge-count prior. Level \\(l \ge 1\\) adds the partition
prior of its \\(B\_{l-1}\\) nodes and the dense multigraph likelihood
\\(\sum\_{r \le s} \log \left(\\!\binom{n_r n_s}{e\_{rs}}\\!\right)\\)
of the group-pair counts. The top level adds the edge-count prior.

`minimize_nested_blockmodel_dl()` follows graph-tool's algorithm. It
starts with one group at every level and \\(\lceil \log_2 N \rceil +
1\\) levels. It then runs the multilevel search of
[`minimize_blockmodel_dl`](https://diogro.github.io/rgraphtool/reference/minimize_blockmodel_dl.md)
on one level at a time, starting at level 0. Every move accounts for the
change in description length at all higher levels. A result is kept only
if it lowers the total description length. After an improvement the
search moves one level up; otherwise it moves one level down, stopping
when it falls below level 0. Refinement then sweeps all levels until the
relative change is at most `epsilon`. The differences from graph-tool
listed in
[`minimize_blockmodel_dl`](https://diogro.github.io/rgraphtool/reference/minimize_blockmodel_dl.md)
apply here too. In addition, if the single-group start violates `B_min`
or `B_max`, the first level-0 result is accepted even when it does not
lower the description length. graph-tool would keep the single group in
that case.

`nested_level_entropy()` returns each level's contribution; they sum to
`nested_description_length()`. `nested_project()` maps vertices to their
groups at a given level. `nested_level_state()` returns that projection
as a flat
[`block_state`](https://diogro.github.io/rgraphtool/reference/block_state.md).
`nested_summary()` lists, for each level up to the first single-node
level, the number of nodes \\(N\\), the number of groups \\(B\\), and
its contribution to the description length. Printing a state shows the
same summary as graph-tool's `print_summary()`.

Nested states keep a reference to the graph and fail with an error if
the graph gains edges or vertices afterwards.

## Value

`minimize_nested_blockmodel_dl()` and `nested_block_state()` return a
nested block state: a list with `levels` (1-based contiguous labels),
`graph`, `degree_corrected`, `degree_prior` and, for fitted states,
`fit`. `fit` holds `description_length` and `trace`, the total
description length after each level search and each refinement pass.

## References

Peixoto, T. P. (2014). Hierarchical block structures and high-resolution
model selection in large networks. *Physical Review X* 4, 011047.
doi:10.1103/PhysRevX.4.011047

## See also

[`plot_hierarchy`](https://diogro.github.io/rgraphtool/reference/plot_hierarchy.md),
[`minimize_blockmodel_dl`](https://diogro.github.io/rgraphtool/reference/minimize_blockmodel_dl.md),
`demo("celegans")`

## Examples

``` r
g <- graph_create(12, directed = FALSE)
graph_add_edges(g, do.call(rbind, lapply(list(1:4, 5:8, 9:12),
                                         function(v) t(combn(v, 2)))))
set.seed(1)
state <- minimize_nested_blockmodel_dl(g)
state
#> <nested block state: degree-corrected, 5 levels, description length 37.2556 nats>
#> l: 0, N: 12, B: 1
#> l: 1, N: 1, B: 1
nested_summary(state)
#>   level  N B description_length
#> 1     0 12 1           37.25555
#> 2     1  1 1            0.00000
nested_project(state, 0)
#>  [1] 1 1 1 1 1 1 1 1 1 1 1 1

fixed <- nested_block_state(g, list(rep(1:3, each = 4), c(1, 1, 1)))
nested_description_length(fixed)
#> [1] 45.29285
```
