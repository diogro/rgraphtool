# Assortative planted-partition stochastic block model

Native C++ implementation of the microcanonical Bayesian
planted-partition stochastic block model (Zhang and Peixoto 2020),
matching graph-tool's `PPBlockState`. Constrains community structure to
assortative (or disassortative) block patterns parameterized only by
within-community edges \\E\_{\mathrm{in}}\\ and between-community edges
\\E\_{\mathrm{out}}\\, and provides declarative model comparison against
the general flat and nested SBMs.

## Usage

``` r
pp_block_state(graph, blocks = NULL, uniform = FALSE,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"))
pp_description_length(state, degree_prior = NULL)
pp_assortativity(state)
minimize_pp_blockmodel_dl(graph, uniform = FALSE, degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"), B_min = 1L, B_max = Inf,
  refine = TRUE, epsilon = 1e-3)
compare_assortative_models(graph,
  models = c("pp", "pp_uniform", "sbm_flat", "sbm_nested"),
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  method = c("map", "posterior"),
  evidence_method = c("variational", "importance"),
  nruns = 1L, nsamples = 200L, thin = 5L, burnin = NULL, mle = TRUE,
  B_min = 1L, B_max = Inf, refine = TRUE)
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
  where each community \\r\\ has its own internal edge count
  \\e\_{rr}\\; if `TRUE`, constrain all communities to share a uniform
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

## Details

Unlike the general SBM, which estimates an arbitrary \\B \times B\\
block affinity matrix, the planted-partition model restricts
off-diagonal interactions to a single global rate while modeling
internal edges either per-community (`uniform = FALSE`) or uniformly
(`uniform = TRUE`). This gives \\O(B)\\ state memory and zero-allocation
incremental moves while remaining free of modularity's resolution limit
and overfitting on random graphs.

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
`degree_corrected`, `degree_prior`, `B`, and optional `fit`.

`pp_description_length()` returns the exact microcanonical description
length in nats.

`pp_assortativity()` returns a list with `E_in`, `E_out`, `E_total`,
`fraction_in`, `B`, effective group count `Be`, and per-group internal
edge counts `err`.

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
