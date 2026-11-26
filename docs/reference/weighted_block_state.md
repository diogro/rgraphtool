# Flat and nested weighted stochastic block models with edge covariates

Construct, evaluate, and infer flat and hierarchical (nested) weighted
stochastic block models (Peixoto 2018) supporting all five graph-tool
edge covariate distributions (`real-exponential`, `real-normal`,
`discrete-geometric`, `discrete-poisson`, and `discrete-binomial`),
multiple simultaneous covariates, and both microcanonical and Bayesian
conjugate hyperpriors.

## Usage

``` r
weighted_block_state(graph, weights, types, params = "microcanonical",
  blocks = NULL, degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"))
weighted_description_length(state, rec = TRUE, rec_dl = TRUE,
  degree_prior = NULL)
minimize_weighted_blockmodel_dl(graph, weights, types,
  params = "microcanonical", degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"), B_min = 1L, B_max = Inf,
  refine = TRUE, epsilon = 1e-3, nthreads = NULL)
nested_weighted_block_state(graph, weights, types, params = "microcanonical",
  levels = NULL, degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"))
nested_weighted_level_entropy(state, degree_prior = NULL)
nested_weighted_description_length(state, degree_prior = NULL)
minimize_nested_weighted_blockmodel_dl(graph, weights, types,
  params = "microcanonical", degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"), B_min = 1L, B_max = Inf,
  refine = TRUE, epsilon = 1e-3, nthreads = NULL)
```

## Arguments

- graph:

  A graph created by
  [`graph_create`](https://diogro.github.io/rgraphtool/reference/graph.md).

- weights:

  Edge covariate values: either a character vector naming one or more
  edge properties stored on `graph`, a numeric vector of length
  `graph_num_edges(graph)`, an \\(E \times K\\) numeric matrix, or a
  list mixing property names and numeric vectors.

- types:

  Character vector of length \\(K\\) naming the covariate model for each
  column of `weights`: `"real-exponential"`, `"real-normal"`,
  `"discrete-geometric"`, `"discrete-poisson"`, or
  `"discrete-binomial"`.

- params:

  Covariate hyperprior specification: `"microcanonical"` (default,
  Peixoto 2018 maximum-entropy microcanonical priors), `"bayesian"`
  (empirical-Bayes conjugate hyperpriors matching graph-tool), a named
  list of hyperparameters (e.g. `list(alpha = 1, beta = 2)` or
  `list(N = 10, alpha = 1, beta = 1)` or
  `list(m0 = 0, k0 = 1, v0 = 1, nu0 = 3)`), or a list of length \\(K\\)
  specifying one entry per covariate.

- blocks:

  Optional 1-based integer group membership vector of length
  `graph_num_vertices(graph)`.

- levels:

  Optional non-empty list of 1-based integer label vectors for each
  hierarchy level.

- degree_corrected:

  Whether to use the degree-corrected (`TRUE`) or non-degree-corrected
  (`FALSE`) topology model.

- degree_prior:

  Degree sequence hyperprior (`"distributed"` or `"uniform"`).

- state:

  A `rgraphtool_weighted_block_state` or
  `rgraphtool_nested_weighted_block_state` object.

- rec:

  Whether to include the edge-covariate log-likelihood term.

- rec_dl:

  Whether to include the global covariate hyperprior description length
  term.

- B_min, B_max:

  Minimum and maximum number of occupied level-0 groups.

- refine:

  Whether to run greedy local refinement sweeps.

- epsilon:

  Convergence tolerance in nats.

- nthreads:

  Optional thread count override (defaults to
  [`rgraphtool_num_threads()`](https://diogro.github.io/rgraphtool/reference/rgraphtool_threads.md)).

## Details

All description lengths match graph-tool 3.9's `WeightedBlockState` and
`NestedBlockState(base_state=WeightedBlockState)` to machine precision
(including the exact continuous quantization scale \\(\epsilon = 0.1
\min\_{w_e \neq 0} \|w_e\|\\) and the nested `real-normal` variance
hyperprior).

Weighted states integrate directly with
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
`rgraphtool_weighted_block_state`; nested constructors and minimizers
return a `rgraphtool_nested_weighted_block_state`. Description length
functions return a numeric scalar (or per-level numeric vector for
`nested_weighted_level_entropy()`).

## References

Peixoto, T. P. (2018). Nonparametric weighted stochastic block models.
*Physical Review E* 97, 012306.

## See also

[`compare_weight_models`](https://diogro.github.io/rgraphtool/reference/compare_weight_models.md),
[`transform_weights`](https://diogro.github.io/rgraphtool/reference/graph.md),
[`graph_set_edge_property`](https://diogro.github.io/rgraphtool/reference/graph.md),
`demo("weighted_sbm")`

## Examples

``` r
set.seed(42)
g <- graph_create(12, directed = FALSE)
edges <- rbind(t(utils::combn(1:6, 2)), t(utils::combn(7:12, 2)), c(6, 7))
w <- c(stats::rnorm(15, 5, 0.5), stats::rnorm(15, -5, 0.5), 0)
graph_add_edges(g, edges, properties = list(w = w))
st <- minimize_weighted_blockmodel_dl(g, weights = "w", types = "real-normal")
st
#> <weighted block state: degree-corrected, real-normal, B = 2, description length 76.6180 nats>
nst <- minimize_nested_weighted_blockmodel_dl(g, weights = "w",
                                              types = "real-normal")
nst
#> <nested weighted block state: degree-corrected, real-normal, 2 levels, B = 2-1, description length 86.0293 nats>
```
