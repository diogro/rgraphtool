# Assortative planted-partition stochastic block model

Native C++ implementation of the microcanonical Bayesian
planted-partition stochastic block model (Zhang and Peixoto 2020),
matching graph-tool's `PPBlockState`. Constrains community structure to
assortative (or disassortative) block patterns parameterized only by
within-community edges \\(E\_{\mathrm{in}}\\) and between-community
edges \\(E\_{\mathrm{out}}\\), optionally paired with exponential-family
or microcanonical edge-weight covariate distributions, and provides
declarative model comparison against the general flat and nested SBMs.

## Usage

``` r
pp_block_state(graph, blocks = NULL, uniform = FALSE,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  weights = NULL, types = NULL, params = "microcanonical",
  layers = NULL, independent = FALSE)
pp_description_length(state, degree_prior = NULL)
pp_assortativity(state)
minimize_pp_blockmodel_dl(graph, uniform = FALSE, degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"), B_min = 1L, B_max = Inf,
  refine = TRUE, epsilon = 1e-3,
  weights = NULL, types = NULL, params = "microcanonical",
  layers = NULL, independent = FALSE)
compare_assortative_models(graph,
  models = c("pp", "pp_uniform", "sbm_flat", "sbm_nested"),
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  method = c("map", "posterior"),
  evidence_method = c("variational", "importance"),
  nruns = 1L, nsamples = 200L, thin = 5L, burnin = NULL, mle = TRUE,
  B_min = 1L, B_max = Inf, refine = TRUE,
  weights = NULL, types = NULL, params = "microcanonical")
```

## Arguments

- graph:

  A graph created by
  [`graph_create`](https://diogro.github.io/rgraphtool/reference/graph.md).

- blocks:

  Optional integer vector of 1-based group memberships of length
  `graph_num_vertices(graph)`. Defaults to a single group.

- uniform:

  If `FALSE` (the default), use the non-uniform planted-partition model
  where each community \\(r\\) has its own internal edge count
  \\(e\_{rr}\\); if `TRUE`, constrain all communities to share a uniform
  internal rate.

- degree_corrected:

  Whether to fit the degree-corrected (`TRUE`) or non-degree-corrected
  (`FALSE`) planted-partition model.

- degree_prior:

  Degree sequence hyperprior (`"distributed"` or `"uniform"`).

- state:

  A `rgraphtool_pp_block_state` created by `pp_block_state()` or
  `minimize_pp_blockmodel_dl()`.

- B_min, B_max:

  Minimum and maximum number of occupied groups.

- refine:

  Whether to run greedy refinement sweeps after merge-split
  minimization.

- epsilon:

  Stopping tolerance in nats for greedy refinement.

- models:

  Character vector naming a subset of `"pp"`, `"pp_uniform"`,
  `"sbm_flat"`, and `"sbm_nested"` to compare.

- method:

  Compare models by `"map"` minimum description length or full Bayesian
  `"posterior"` marginal evidence.

- evidence_method:

  Evidence estimator forwarded to
  [`model_evidence`](https://diogro.github.io/rgraphtool/reference/compare_model_classes.md)
  when `method = "posterior"`.

- nruns:

  Number of independent minimization restarts per model.

- nsamples, thin, burnin, mle:

  Posterior sampling and entropy estimation controls forwarded to
  [`model_evidence`](https://diogro.github.io/rgraphtool/reference/compare_model_classes.md).

- weights:

  Optional edge covariate specification (edge property name(s), vector
  of length `graph_num_edges(graph)`, or matrix/data.frame with one
  column per covariate) using the same distributions and conventions as
  [`weighted_block_state`](https://diogro.github.io/rgraphtool/reference/weighted_block_state.md).
  When provided, `compare_assortative_models()` compares weighted
  planted-partition models against weighted flat and nested SBMs
  ([`minimize_weighted_blockmodel_dl`](https://diogro.github.io/rgraphtool/reference/weighted_block_state.md)
  and
  [`minimize_nested_weighted_blockmodel_dl`](https://diogro.github.io/rgraphtool/reference/weighted_block_state.md)).

- types:

  Optional character vector of covariate families corresponding to
  `weights` (`"real-normal"`, `"real-exponential"`,
  `"discrete-poisson"`, `"discrete-geometric"`, or
  `"discrete-binomial"`). Required when `weights` is supplied.

- params:

  Hyperprior specification for `weights`: `"microcanonical"` (default)
  or explicit Bayesian prior hyperparameters as in
  [`weighted_block_state`](https://diogro.github.io/rgraphtool/reference/weighted_block_state.md).

- layers:

  Optional discrete edge layer labels: an edge property name or a vector
  of length `graph_num_edges(graph)` (as in
  [`layered_block_state`](https://diogro.github.io/rgraphtool/reference/layered_block_state.md)).
  Enables flat multi-layer planted-partition inference, optionally
  combined with `weights`.

- independent:

  Used when `layers` is provided: if `FALSE` (default), fit the discrete
  edge-covariate PP model with pooled degrees and pooled covariate
  buckets across layers; if `TRUE`, fit conditionally independent PP
  layers with layer-specific degrees, active group sizes, and (when
  weighted) per-layer covariate distributions.

## Details

Unlike the general SBM, which estimates an arbitrary \\(B \times B\\)
block affinity matrix, the planted-partition model restricts
off-diagonal interactions to a single global rate while modeling
internal edges either per-community (`uniform = FALSE`) or uniformly
(`uniform = TRUE`). This gives \\(O(B)\\) state memory and
zero-allocation incremental moves while remaining free of modularity's
resolution limit and overfitting on random graphs.

**Weighted planted-partition models and caveats:** When `weights` and
`types` are supplied, each active edge bucket receives its own weight
distribution under the exact conjugate/microcanonical likelihoods of
[`weighted_block_state`](https://diogro.github.io/rgraphtool/reference/weighted_block_state.md)
(2 buckets `in`/`out` when `uniform = TRUE`; \\(B + 1\\) buckets
\\(\\{rr\\}\\) and `out` when `uniform = FALSE`):

- *Pooled inter-module background:* All between-community edge weights
  are pooled into a single `out` distribution rather than modeling
  pair-specific \\((r, s)\\) block weight differences. If specific pairs
  of modules have distinct positive vs. negative weight patterns (e.g.,
  inter-module antagonism in a gene regulatory network), use
  `compare_assortative_models(..., weights = ..., types = ...)` to test
  whether `"sbm_flat"` or `"sbm_nested"` has higher Bayesian evidence.

- *Mean-ordering / sign symmetry (especially at \\(B = 2\\)):* With
  signed covariates such as `"real-normal"`, the covariate likelihood
  tests whether the within-community and between-community weight
  distributions *differ*, without enforcing \\(\mu\_{\mathrm{in}} \>
  \mu\_{\mathrm{out}}\\). On strongly anti-correlated bipartitions
  (particularly at \\(B = 2\\)), inspect `pp_assortativity()`
  (`fraction_in`, `weight_in_mean`, `weight_out_mean`) to confirm
  whether a recovered partition is assortative or disassortative.

- *Extension beyond upstream graph-tool 3.9:* Upstream `graph-tool` 3.9
  `PPBlockState` does not support edge covariates (`recs`) or edge
  layers; `rgraphtool` extends PP to both (flat, single-level) layers
  and covariates using the exact `graph-tool` microcanonical/conjugate
  likelihoods. At \\(B = 1\\), the weighted/layered PP description
  length matches
  [`weighted_block_state`](https://diogro.github.io/rgraphtool/reference/weighted_block_state.md)
  or
  [`layered_block_state`](https://diogro.github.io/rgraphtool/reference/layered_block_state.md)
  identically, and
  [`compare_layered_models()`](https://diogro.github.io/rgraphtool/reference/compare_layered_models.md)
  supports direct comparison against general layered SBMs.

Planted-partition states work directly with
[`block_membership`](https://diogro.github.io/rgraphtool/reference/block_state.md),
[`block_sizes`](https://diogro.github.io/rgraphtool/reference/block_state.md),
[`mcmc_sweep`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md),
[`mcmc_equilibrate`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md),
[`mcmc_anneal`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md),
[`sample_posterior`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md),
and
[`model_evidence`](https://diogro.github.io/rgraphtool/reference/compare_model_classes.md).

## Value

`pp_block_state()` and `minimize_pp_blockmodel_dl()` return a
`rgraphtool_pp_block_state` list with `blocks`, `graph`, `uniform`,
`degree_corrected`, `degree_prior`, `B`, and optional `fit` (plus
`weights`, `types`, `params`, `weighted`, `layers`, `layer_levels`, `C`,
`independent`, and `layered`).

`pp_description_length()` returns the exact microcanonical description
length in nats.

`pp_assortativity()` returns a list with `E_in`, `E_out`, `E_total`,
`fraction_in`, `B`, effective group count `Be`, and per-group internal
edge counts `err` (plus `weight_in_mean` and `weight_out_mean` when
`state` is weighted, and per-layer fields `E_in_layer`, `E_out_layer`,
`E_total_layer`, `fraction_in_layer`, `err_layer`,
`weight_in_mean_layer`, and `weight_out_mean_layer` when `state` is
layered).

`compare_assortative_models()` returns a
`rgraphtool_assortative_model_comparison` list with ranked `summary`,
`best_model`, `best_state`, `states`, optional `evidences`, and
`log_bayes_factors`.

## References

Zhang, L. and Peixoto, T. P. (2020). Statistical inference of
assortative community structures. *Physical Review Research* 2, 043271.
doi:10.1103/PhysRevResearch.2.043271

## See also

[`minimize_blockmodel_dl`](https://diogro.github.io/rgraphtool/reference/minimize_blockmodel_dl.md),
[`minimize_nested_blockmodel_dl`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md),
[`mcmc_sweep`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md),
[`model_evidence`](https://diogro.github.io/rgraphtool/reference/compare_model_classes.md),
`demo("assortative_sbm")`

## Examples

``` r
g <- graph_create(12, directed = FALSE)
graph_add_edges(g, rbind(t(combn(1:6, 2)), t(combn(7:12, 2)), c(1, 7)))
set.seed(1)
st <- minimize_pp_blockmodel_dl(g)
st
#> <planted-partition block state: degree-corrected, non-uniform assortativity, B = 2 (Be = 2.00), E_in/E = 96.8%, DL = 55.1833 nats>
pp_assortativity(st)
#> $E_in
#> [1] 30
#> 
#> $E_out
#> [1] 1
#> 
#> $E_total
#> [1] 31
#> 
#> $fraction_in
#> [1] 0.9677419
#> 
#> $B
#> [1] 2
#> 
#> $Be
#> [1] 2
#> 
#> $err
#> [1] 15 15
#> 
```
