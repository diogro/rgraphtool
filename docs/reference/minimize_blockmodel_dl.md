# Fit a flat stochastic block model by minimum description length

`minimize_blockmodel_dl()` is the R counterpart of graph-tool's
`minimize_blockmodel_dl()` for the flat `BlockState`. It runs a native
multilevel agglomerative search followed by refinement, and returns the
partition with the smallest description length it found.

## Usage

``` r
minimize_blockmodel_dl(graph, degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"), B_min = 1, B_max = Inf,
  refine = TRUE, epsilon = 1e-3, init_blocks = NULL, nthreads = NULL,
  control = list())
refine_blockmodel_dl(graph, state, degree_prior = NULL, epsilon = 1e-3,
  nthreads = NULL, control = list())
compare_blockmodels(graph, degree_prior = c("distributed", "uniform"),
  nested = FALSE, nruns = 1L, ...)
```

## Arguments

- graph:

  A graph created by
  [`graph_create()`](https://diogro.github.io/rgraphtool/reference/graph.md).

- degree_corrected:

  Whether to fit the degree-corrected model.

- degree_prior:

  Degree-sequence prior, `"distributed"` (graph-tool's default) or
  `"uniform"`. For `refine_blockmodel_dl()`, `NULL` uses the prior
  stored in `state`.

- B_min, B_max:

  Bounds on the number of nonempty groups.

- refine:

  Whether to run the refinement loop after the multilevel search.

- epsilon:

  Refinement stops when a pass changes the description length by a
  relative amount of at most `epsilon`.

- init_blocks:

  Optional warm-start partition (either a 1-based integer vector of
  length `graph_num_vertices(graph)` or an existing
  `rgraphtool_block_state`), skipping initialization from singleton
  groups.

- nthreads:

  Optional positive integer number of parallel worker threads (or
  `"auto"`). `NULL` uses
  [`rgraphtool_num_threads`](https://diogro.github.io/rgraphtool/reference/rgraphtool_threads.md).

- control:

  Named list of search settings; see Details.

- state:

  A state created by
  [`block_state()`](https://diogro.github.io/rgraphtool/reference/block_state.md)
  or a fitting function.

- nested:

  Whether `compare_blockmodels()` fits nested models with
  [`minimize_nested_blockmodel_dl`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md)
  instead of flat ones.

- nruns:

  Number of fits per model; the one with the smallest description length
  is kept.

- ...:

  Passed to the fitting function.

## Details

The search follows graph-tool 3.9's `multilevel_mcmc_sweep()` at zero
temperature (\\(\beta = \infty\\)) over all groups:

1.  Every vertex starts in its own group. Single-vertex sweeps then move
    vertices into neighbouring groups, starting at \\(\beta =\\)
    `init_beta` for `init_min_iter_beta` sweeps and continuing at zero
    temperature for at least `init_min_iter` sweeps, until a sweep
    reduces the number of groups by less than a factor of `init_r`.

2.  A bisection search with random midpoints looks for the number of
    groups \\(B\\) that minimizes the description length. To reach a
    smaller \\(B\\) from a cached partition, the search alternates
    merges (each group tries `merge_sweeps` neighbour-proposed merge
    partners, and the best merges are applied until \\(B\\) falls by a
    factor of `r`) with up to `mh_sweeps` zero-temperature single-vertex
    sweeps. Every partition visited is cached by \\(B\\).

3.  Refinement repeatedly merges groups when that lowers the description
    length and sweeps vertices among the existing groups. It stops after
    a pass whose relative change is at most `epsilon`.

Moves never increase the description length after the initial
`init_beta` sweeps. The result is a local optimum, not a guaranteed
global one, so run several seeds and keep the smallest description
length, as the graph-tool documentation recommends.

`control` accepts `r` (0.9), `merge_sweeps` (10), `mh_sweeps` (10),
`init_r` (0.99), `init_min_iter` (8), `init_min_iter_beta` (4),
`init_beta` (1), `c` (0.5, the random proposal weight used when a group
has no neighbour-proposed merge partner), `max_refine` (1000),
`bisection` (`TRUE`; set to `FALSE` for strictly agglomerative geometric
coarsening without bisection backtracking on large networks), `nthreads`
(overrides
[`rgraphtool_num_threads`](https://diogro.github.io/rgraphtool/reference/rgraphtool_threads.md)
when `nthreads = NULL`), and `verbose` (`FALSE`).

When `nthreads > 1`, candidate merge evaluations in `merge_sweep()` and
zero-temperature vertex proposal screening in `mh_sweep()` run in
parallel across a portable C++23 thread pool, while all state mutations
and priority-queue merge commits execute deterministically on the main
thread. For very large networks (\\(N \ge 10^4\\)), combining `nthreads`
with strictly agglomerative search
(`control = list(bisection = FALSE, r = 0.5)`) keeps only a single
snapshot in memory and reduces runtime by an order of magnitude.

`compare_blockmodels()` fits the non-degree-corrected (NDC) and
degree-corrected (DC) models and reports which has the smaller
description length \\(\Sigma\\). Both description lengths include all
priors, so with equal prior probabilities for the two models the
posterior odds ratio is \$\$\Lambda = \frac{P(b, H\_{NDC} \mid A)}{P(b,
H\_{DC} \mid A)} = e^{-(\Sigma\_{NDC} - \Sigma\_{DC})},\$\$ as in the
graph-tool cookbook. `log_posterior_odds` is \\(\ln \Lambda =
\Sigma\_{DC} - \Sigma\_{NDC}\\), so negative values favour the
degree-corrected model. (The cookbook prints this difference with the
opposite sign for the football network.) Because the fits are stochastic
local optima, increase `nruns` before relying on small differences.

## Value

`minimize_blockmodel_dl()` and `refine_blockmodel_dl()` return a block
state (see
[`block_state`](https://diogro.github.io/rgraphtool/reference/block_state.md))
with labels `1, ..., B`. Its `fit` element holds `description_length`.
For `minimize_blockmodel_dl()` it also holds
`multilevel_description_length`, the score before refinement, and
`trace`, a data frame of every cached \\((B, \mathrm{DL})\\) pair in the
order visited. The description length of `refine_blockmodel_dl()`'s
result is never higher than that of the input state.

`compare_blockmodels()` returns a list with `comparison`, a data frame
of models and description lengths; `log_posterior_odds`; `best_model`;
`best_state`; `states`, both fitted states; and `nested`. Printing it
shows the description lengths and \\(\ln \Lambda\\).

## References

Peixoto, T. P. (2014). Efficient Monte Carlo and greedy heuristic for
the inference of stochastic block models. *Physical Review E* 89,
012804. doi:10.1103/PhysRevE.89.012804

Peixoto, T. P. (2017). Nonparametric Bayesian inference of the
microcanonical stochastic block model. *Physical Review E* 95, 012317.
doi:10.1103/PhysRevE.95.012317

## See also

[`block_state`](https://diogro.github.io/rgraphtool/reference/block_state.md)

## Examples

``` r
g <- graph_create(12, directed = FALSE)
graph_add_edges(g, do.call(rbind, lapply(list(1:6, 7:12),
                                         function(v) t(combn(v, 2)))))
set.seed(1)
fit <- minimize_blockmodel_dl(g)
block_membership(fit)
#>  [1] 1 1 1 1 1 1 2 2 2 2 2 2
block_description_length(fit)
#> [1] 47.38183
head(fit$fit$trace)
#>    B  entropy
#> 1  1 50.94809
#> 2 12 83.49776
#> 3  2 47.38183
compare_blockmodels(g)$comparison
#>                  model description_length
#> 1 non-degree-corrected           48.83070
#> 2     degree-corrected           47.38183
```
