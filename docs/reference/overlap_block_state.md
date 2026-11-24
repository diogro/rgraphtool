# Flat and hierarchical overlapping stochastic block models

Native C++ implementation of Peixoto's (2015) microcanonical overlapping
stochastic block model (`OverlapBlockState`) and its hierarchical
extension (`NestedBlockState` with an overlapping base level). Each edge
endpoint (half-edge) carries its own group label, allowing vertices to
participate in multiple communities simultaneously.

## Usage

``` r
overlap_block_state(graph, half_edges = NULL, degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"), gen_overlap = TRUE)
overlap_description_length(state, degree_prior = NULL)
overlap_groups(state)
overlap_majority_groups(state)
minimize_overlap_blockmodel_dl(graph, degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"), gen_overlap = TRUE,
  B_min = 1L, B_max = Inf, refine = TRUE, epsilon = 1e-3)
nested_overlap_block_state(graph, half_edges, upper_levels,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  gen_overlap = TRUE)
nested_overlap_level_entropy(state, degree_prior = NULL)
nested_overlap_description_length(state, degree_prior = NULL)
minimize_nested_overlap_blockmodel_dl(graph, degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"), gen_overlap = TRUE,
  B_min = 1L, B_max = Inf, refine = TRUE, epsilon = 1e-3)
overlap_mcmc_sweep(state, beta = 1, niter = 1L, c = 0.5, d = 0.01,
  thin = 1L, collect = FALSE, B_min = 1L, B_max = Inf)
```

## Arguments

- graph:

  A graph created by
  [`graph_create`](https://diogro.github.io/rgraphtool/reference/graph.md).

- half_edges:

  An `E x 2` integer matrix of positive group labels for the source and
  target half-edges of each row in
  [`graph_edges`](https://diogro.github.io/rgraphtool/reference/graph.md)`(graph)`,
  or `NULL` to place all half-edges in a single group.

- upper_levels:

  Optional list of 1-based integer label vectors defining hierarchical
  partitions above the overlapping base level.

- degree_corrected:

  Logical indicating whether to use the degree-corrected overlapping
  SBM.

- degree_prior:

  Degree sequence hyperprior: `"distributed"` (default, Peixoto's
  hyperprior) or `"uniform"`.

- gen_overlap:

  Logical indicating whether to use the generative overlapping mixture
  and degree hyperpriors (default `TRUE`, matching graph-tool 3.9).

- state:

  An object returned by `overlap_block_state()`,
  `minimize_overlap_blockmodel_dl()`, `nested_overlap_block_state()`, or
  `minimize_nested_overlap_blockmodel_dl()`.

- B_min, B_max:

  Minimum and maximum number of occupied level-0 groups.

- refine:

  Whether to run greedy half-edge, vertex, edge-coupled, merge, split,
  and bridge-dissolution refinement sweeps.

- epsilon:

  Relative entropy convergence tolerance for refinement.

- beta:

  Inverse temperature for Metropolis-Hastings sweeps (`1` for exact
  Bayesian posterior sampling; `Inf` for greedy moves).

- niter:

  Number of sweeps to perform.

- c, d:

  Proposal mixing parameter and empty-group proposal probability.

- thin:

  Keep every `thin`-th sweep when `collect = TRUE`.

- collect:

  Whether to collect vertex majority-group partitions and description
  lengths during MCMC sweeps.

## Details

Exact flat and nested overlapping description lengths match graph-tool
3.9 to machine precision for directed and undirected graphs, with or
without parallel edges, self-loops, isolated vertices, degree
correction, and distributed or uniform degree priors.

Minimization combines non-overlapping seeding, vertex- and edge-level
agglomeration, sequential bridge-group dissolution, and greedy
half-edge, vertex, and edge-coupled refinement sweeps. Per-vertex mixed
memberships returned by `overlap_groups(state)$matrix` can be plotted
directly as pie-chart nodes via
[`plot_marginals`](https://diogro.github.io/rgraphtool/reference/plot_marginals.md).

## Value

`overlap_block_state()` and `minimize_overlap_blockmodel_dl()` return a
`rgraphtool_overlap_block_state` object containing `half_edges`,
`overlap`, and model settings. `nested_overlap_block_state()` and
`minimize_nested_overlap_blockmodel_dl()` return a
`rgraphtool_nested_overlap_block_state` object adding `upper_levels`.
`overlap_groups()` returns a list with per-vertex `groups`,
`out_degrees`, `in_degrees`, `total_degrees`, half-edge count `matrix`,
`majority` group vector, and occupied group count `B`.
`overlap_mcmc_sweep()` returns a list with updated `state`, `entropy`,
`dS`, `nattempts`, `nmoves`, and (when `collect = TRUE`)
`majority_samples` and `entropies`.

## References

Peixoto, T. P. (2015). Model selection and hypothesis testing for
large-scale network models with overlapping groups. *Physical Review X*
5, 011033. doi:10.1103/PhysRevX.5.011033

## See also

[`compare_model_classes`](https://diogro.github.io/rgraphtool/reference/compare_model_classes.md),
[`model_evidence`](https://diogro.github.io/rgraphtool/reference/compare_model_classes.md),
[`plot_marginals`](https://diogro.github.io/rgraphtool/reference/plot_marginals.md),
`demo("model_class_selection")`

## Examples

``` r
g <- graph_create(n = 6L, directed = FALSE)
graph_add_edges(g, matrix(c(1, 2, 1, 3, 2, 3, 3, 4, 3, 5, 4, 5),
                          ncol = 2L, byrow = TRUE))
he <- matrix(c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2), ncol = 2L, byrow = TRUE)
st <- overlap_block_state(g, half_edges = he)
overlap_description_length(st)
#> [1] 16.66993
overlap_groups(st)$matrix
#>      [,1] [,2]
#> [1,]    2    0
#> [2,]    2    0
#> [3,]    2    2
#> [4,]    0    2
#> [5,]    0    2
#> [6,]    0    0
```
