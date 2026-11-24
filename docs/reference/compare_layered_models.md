# Compare edge-covariate and independent-layers SBM variants

Declarative, callback-free model selection between the discrete
edge-covariate layered SBM (`independent = FALSE`) and the conditionally
independent layers SBM (`independent = TRUE`), optionally across
degree-corrected and non-degree-corrected variants and using either MAP
description lengths or full Bayesian posterior marginal evidence.

## Usage

``` r
compare_layered_models(graph, layers,
  models = c("edge_covariates", "independent_layers"), nested = TRUE,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  method = c("map", "posterior"),
  evidence_method = c("variational", "importance"), nruns = 1L,
  nsamples = 200L, thin = 5L, burnin = NULL, mle = TRUE, B_min = 1L,
  B_max = Inf, refine = TRUE)
```

## Arguments

- graph:

  A graph created by
  [`graph_create`](https://diogro.github.io/rgraphtool/reference/graph.md).

- layers:

  Discrete edge layer labels: an edge property name or a vector of
  length `graph_num_edges(graph)`.

- models:

  Character vector naming a subset of `"edge_covariates"`,
  `"independent_layers"`, `"edge_covariates_dc"`,
  `"edge_covariates_ndc"`, `"independent_layers_dc"`, and
  `"independent_layers_ndc"`. The unqualified names use
  `degree_corrected`.

- nested:

  Whether to fit nested (`TRUE`, default) or flat (`FALSE`) layered
  SBMs.

- degree_corrected:

  Degree correction setting used for `"edge_covariates"` and
  `"independent_layers"`.

- degree_prior:

  Degree sequence hyperprior (`"distributed"` or `"uniform"`).

- method:

  Rank candidate models by `"map"` minimum description length or full
  Bayesian `"posterior"` marginal evidence.

- evidence_method:

  Evidence estimator forwarded to
  [`model_evidence`](https://diogro.github.io/rgraphtool/reference/compare_model_classes.md)
  when `method = "posterior"`.

- nruns:

  Number of independent minimization restarts per candidate.

- nsamples, thin, burnin, mle:

  Posterior sampling arguments used when `method = "posterior"`.

- B_min, B_max:

  Bounds on the number of level-0 groups.

- refine:

  Whether to run local refinement sweeps during minimization.

## Value

A `rgraphtool_layered_model_comparison` list containing:

- `summary`: ranked data frame with `model`, `independent`,
  `degree_corrected`, `B`, `description_length`, `posterior_entropy`,
  `neg_log_evidence`, `delta`, and `posterior_weight`.

- `best_model`: name of the winning layered model.

- `best_state`: fitted layered state for the winning model.

- `states`: named list of fitted states ordered from best to worst.

- `evidences`: named list of `rgraphtool_model_evidence` objects (when
  `method = "posterior"`).

- `log_bayes_factors`: pairwise matrix of log Bayes factors.

## References

Peixoto, T. P. (2015). Inferring the mesoscale structure of layered,
edge-valued, and time-varying networks. *Physical Review E* 92, 042807.

## See also

[`layered_block_state`](https://diogro.github.io/rgraphtool/reference/layered_block_state.md),
[`model_evidence`](https://diogro.github.io/rgraphtool/reference/compare_model_classes.md),
`demo("layered_sbm")`

## Examples

``` r
set.seed(42)
g <- graph_create(12, directed = FALSE)
edges <- rbind(t(utils::combn(1:6, 2)), t(utils::combn(7:12, 2)))
lay <- c(rep("layer_a", 15), rep("layer_b", 15))
graph_add_edges(g, edges)
cmp <- compare_layered_models(g, layers = lay, nested = FALSE)
cmp
#> <layered SBM model comparison (flat, 2 layers, MAP)>
#>               model B description_length  delta posterior_weight
#>  independent_layers 1             36.415  0.000          1.0e+00
#>     edge_covariates 2             54.368 17.953          1.6e-08
#> Best layered model: independent_layers
```
