# Flat and nested layered stochastic block models

Construct, evaluate, and infer flat and hierarchical (nested) layered
stochastic block models (Peixoto 2015), supporting both the
conditionally independent layers model (`independent = TRUE`) and the
discrete edge-covariate model (`independent = FALSE`).

## Usage

``` r
layered_block_state(graph, layers, independent = TRUE, blocks = NULL,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"))
layered_description_length(state, lmemb = TRUE, degree_prior = NULL)
minimize_layered_blockmodel_dl(graph, layers, independent = TRUE,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  B_min = 1L, B_max = Inf, refine = TRUE, epsilon = 1e-3)
nested_layered_block_state(graph, layers, independent = TRUE, levels = NULL,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"))
nested_layered_level_entropy(state, degree_prior = NULL)
nested_layered_description_length(state, degree_prior = NULL)
minimize_nested_layered_blockmodel_dl(graph, layers, independent = TRUE,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  B_min = 1L, B_max = Inf, refine = TRUE, epsilon = 1e-3)
```

## Arguments

- graph:

  A graph created by
  [`graph_create`](https://diogro.github.io/rgraphtool/reference/graph.md).

- layers:

  Discrete edge layer labels: either a character scalar naming an edge
  property stored on `graph` or a numeric, integer, character, factor,
  or logical vector of length `graph_num_edges(graph)`.

- independent:

  If `TRUE` (default), fit the conditionally independent layers model
  (each layer has its own active vertices, degree sequence, and block
  edge counts conditioned on the shared partition). If `FALSE`, fit the
  edge-covariate model (a single collapsed topology and degree sequence
  with categorical layer labels per block pair).

- blocks:

  Optional 1-based integer group membership vector of length
  `graph_num_vertices(graph)`.

- levels:

  Optional non-empty list of 1-based integer label vectors for each
  hierarchy level.

- degree_corrected:

  Whether to use degree correction (`TRUE`) or the non-degree-corrected
  model (`FALSE`).

- degree_prior:

  Degree sequence hyperprior (`"distributed"` or `"uniform"`).

- state:

  A `rgraphtool_layered_block_state` or
  `rgraphtool_nested_layered_block_state` object.

- lmemb:

  Whether to include the layer-membership description length term
  (applicable when `independent = TRUE`).

- B_min, B_max:

  Minimum and maximum number of occupied level-0 groups.

- refine:

  Whether to run local refinement sweeps during minimization.

- epsilon:

  Convergence tolerance in nats.

## Details

All description lengths match graph-tool 3.9's `LayeredBlockState` and
`NestedBlockState(base_state=LayeredBlockState)` to machine precision
for directed and undirected graphs (including parallel edges within and
across layers, self-loops, isolated vertices, and blocks absent from
individual layers).

Layered states integrate directly with
[`block_membership`](https://diogro.github.io/rgraphtool/reference/block_state.md),
[`block_sizes`](https://diogro.github.io/rgraphtool/reference/block_state.md),
[`nested_levels`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md),
[`nested_project`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md),
[`nested_summary`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md),
[`plot_hierarchy`](https://diogro.github.io/rgraphtool/reference/plot_hierarchy.md),
[`mcmc_sweep`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md),
[`mcmc_equilibrate`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md),
[`mcmc_anneal`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md),
[`sample_posterior`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md),
and
[`model_evidence`](https://diogro.github.io/rgraphtool/reference/compare_model_classes.md).

## Value

Flat constructors and minimizers return a
`rgraphtool_layered_block_state`; nested constructors and minimizers
return a `rgraphtool_nested_layered_block_state`. Description length
functions return a numeric scalar (or per-level numeric vector for
`nested_layered_level_entropy()`).

## References

Peixoto, T. P. (2015). Inferring the mesoscale structure of layered,
edge-valued, and time-varying networks. *Physical Review E* 92, 042807.

## See also

[`compare_layered_models`](https://diogro.github.io/rgraphtool/reference/compare_layered_models.md),
[`model_evidence`](https://diogro.github.io/rgraphtool/reference/compare_model_classes.md),
[`graph_set_edge_property`](https://diogro.github.io/rgraphtool/reference/graph.md),
`demo("layered_sbm")`

## Examples

``` r
set.seed(42)
g <- graph_create(10, directed = FALSE)
edges <- rbind(t(utils::combn(1:5, 2)), t(utils::combn(6:10, 2)),
               c(1, 6), c(2, 7), c(3, 8))
lay <- c(rep("within", 20), rep("across", 3))
graph_add_edges(g, edges)
st <- minimize_layered_blockmodel_dl(g, layers = lay, independent = TRUE)
st
#> <layered block state: degree-corrected, 2 independent layers, B = 1, description length 48.3846 nats>
nst <- minimize_nested_layered_blockmodel_dl(g, layers = lay,
                                             independent = TRUE)
nst
#> <nested layered block state: degree-corrected, 2 independent layers, 2 levels, B = 1-1, description length 48.3846 nats>
```
