# Posterior model evidence and SBM model class selection

Declarative, callback-free Bayesian model class selection across flat,
hierarchical (nested), degree-corrected, non-degree-corrected, and
overlapping stochastic block models. Computes both MAP description
lengths and full posterior marginal log-evidences \\(\ln P(A \|
\mathcal{M})\\), posterior model weights, and pairwise log Bayes
factors.

## Usage

``` r
model_evidence(state, nsamples = 200L, thin = 5L, burnin = NULL,
  beta = 1, mle = TRUE, method = c("variational", "importance"),
  samples = NULL, entropies = NULL, multiflip = TRUE,
  c = 0.5, d = 0.01, B_min = 1L, B_max = Inf)
compare_model_classes(graph,
  models = c("flat_ndc", "flat_dc", "nested_ndc", "nested_dc",
    "overlap_flat_ndc", "overlap_flat_dc",
    "overlap_nested_ndc", "overlap_nested_dc"),
  degree_prior = c("distributed", "uniform"),
  method = c("map", "posterior"),
  evidence_method = c("variational", "importance"),
  nruns = 1L, nsamples = 200L, thin = 5L, burnin = NULL, mle = TRUE,
  B_min = 1L, B_max = Inf, refine = TRUE)
```

## Arguments

- state:

  A fitted `rgraphtool_block_state`, `rgraphtool_nested_block_state`,
  `rgraphtool_overlap_block_state`,
  `rgraphtool_nested_overlap_block_state`,
  `rgraphtool_weighted_block_state`,
  `rgraphtool_nested_weighted_block_state`,
  `rgraphtool_layered_block_state`,
  `rgraphtool_nested_layered_block_state`, `rgraphtool_pp_block_state`,
  `rgraphtool_ranked_block_state`, or
  `rgraphtool_nested_ranked_block_state` object.

- graph:

  A graph created by
  [`graph_create`](https://diogro.github.io/rgraphtool/reference/graph.md).

- models:

  Character vector naming a subset of the 8 SBM model classes:
  `"flat_ndc"`, `"flat_dc"`, `"nested_ndc"`, `"nested_dc"`,
  `"overlap_flat_ndc"`, `"overlap_flat_dc"`, `"overlap_nested_ndc"`, and
  `"overlap_nested_dc"`.

- degree_prior:

  Degree sequence hyperprior (`"distributed"` or `"uniform"`).

- method:

  For `compare_model_classes()`, rank models by `"map"` minimum
  description length or full Bayesian `"posterior"` marginal evidence.
  For `model_evidence()`, choose between `"variational"` (\\(\ln P(A) =
  -\langle \Sigma \rangle + H(b)\\)) and `"importance"` (Peixoto 2021
  Eq. 29 mode-mixture importance sampling).

- evidence_method:

  Evidence estimator forwarded to `model_evidence()` when
  `method = "posterior"`.

- nruns:

  Number of independent minimization restarts per model class, retaining
  the lowest description length fit.

- nsamples:

  Number of thinned posterior samples to collect.

- thin:

  Number of MCMC sweeps between collected samples.

- burnin:

  Number of burn-in sweeps before sampling, or `NULL` to equilibrate
  automatically (flat/nested) or use 30 sweeps (overlap).

- beta:

  Inverse temperature (default `1`).

- mle:

  Whether to use MLE or Dirichlet(1) posterior mean marginals for
  posterior entropy and mode log-probabilities.

- samples, entropies:

  Optional pre-collected posterior partition samples and matching
  description lengths (e.g. from
  [`sample_posterior`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md)),
  avoiding redundant MCMC runs.

- multiflip:

  Whether to use merge-split MCMC moves for flat/nested states.

- c, d:

  MCMC proposal mixing and empty-group parameters.

- B_min, B_max:

  Minimum and maximum number of level-0 groups.

- refine:

  Whether to run refinement sweeps during minimization.

## Details

For a given model class, `model_evidence()` equilibrates an MCMC chain,
collects posterior partition samples and exact description lengths
\\(\Sigma(b^{(m)})\\) natively in C++, clusters the samples into
posterior modes via
[`mode_cluster_state`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md),
and computes both the variational evidence \\(\ln P(A) = -\langle \Sigma
\rangle + H\_{\mathrm{posterior}}\\) (matching graph-tool's model class
selection cookbook) and Peixoto's (2021) mode-mixture importance
estimate \\(-\ln P(A) = \ln\\!\left(\frac{1}{M}\sum\_{m=1}^M
\exp(\Sigma(b^{(m)}) + \ln q(b^{(m)}))\right)\\).

Note that for overlapping states, the posterior partition entropy is
estimated from per-vertex majority-group samples rather than full
half-edge trajectories, providing a conservative lower bound on the
half-edge posterior entropy.

## Value

`model_evidence()` returns a `rgraphtool_model_evidence` list with
`log_evidence`, `neg_log_evidence`, `variational_log_evidence`,
`importance_log_evidence`, `mean_dl`, `sd_dl`, `map_dl`,
`posterior_entropy`, `cluster`, `samples`, `entropies`, and `state`.

`compare_model_classes()` returns a `rgraphtool_model_class_comparison`
list with ranked `summary` data frame, `best_model`, `best_state`, named
list of fitted `states`, optional `evidences`, and pairwise
`log_bayes_factors` matrix.

## References

Peixoto, T. P. (2015). Model selection and hypothesis testing for
large-scale network models with overlapping groups. *Physical Review X*
5, 011033.

Peixoto, T. P. (2021). Revealing consensus and dissensus between network
partitions. *Physical Review X* 11, 021003.

## See also

[`overlap_block_state`](https://diogro.github.io/rgraphtool/reference/overlap_block_state.md),
[`mode_cluster_state`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md),
[`sample_posterior`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md),
`demo("model_class_selection")`

## Examples

``` r
g <- graph_create(n = 12L, directed = FALSE)
edges <- rbind(
  t(utils::combn(1:6, 2)),
  t(utils::combn(7:12, 2)),
  c(6, 7)
)
graph_add_edges(g, edges)
set.seed(42)
cmp <- compare_model_classes(g, models = c("flat_ndc", "flat_dc",
                                           "nested_ndc", "overlap_flat_ndc"))
cmp
#> <SBM model class comparison (MAP description length, distributed degree prior)>
#>             model B levels description_length delta posterior_weight
#>  overlap_flat_ndc 2      1             50.159 0.000           0.8620
#>          flat_ndc 2      1             52.477 2.318           0.0849
#>        nested_ndc 2      5             53.170 3.011           0.0424
#>           flat_dc 2      1             54.553 4.394           0.0106
#> Best model: overlap_flat_ndc
```
