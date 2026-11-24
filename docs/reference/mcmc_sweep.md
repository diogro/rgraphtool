# MCMC posterior sampling, equilibration, and simulated annealing for SBMs

Native C++ Metropolis-Hastings and merge-split Markov chain Monte Carlo
samplers for flat and nested microcanonical stochastic block models,
targeting the exact posterior \\P(\mathbf{b} \mid A) \propto \exp(-\beta
\Sigma(\mathbf{b}))\\. These are the R counterparts of graph-tool's
`mcmc_sweep()`, `multiflip_mcmc_sweep()`, `mcmc_equilibrate()`, and
`mcmc_anneal()`.

## Usage

``` r
mcmc_sweep(state, beta = 1, niter = 1L, c = NULL, d = 0.01,
  multiflip = TRUE, B_min = 1L, B_max = Inf, psingle = NULL,
  psplit = 1, pmerge = 1, pmergesplit = 1, gibbs_sweeps = 10L,
  ls_shuffle = TRUE)
mcmc_equilibrate(state, wait = 1000L, nbreaks = 2L, max_niter = Inf,
  force_niter = NULL, epsilon = 0, multiflip = TRUE, mcmc_args = list(),
  history = FALSE, callback = NULL, collect = FALSE, thin = 1L,
  verbose = FALSE)
mcmc_anneal(state, beta_range = c(1, 10), niter = 100L,
  history = FALSE, mcmc_equilibrate_args = list(), verbose = FALSE)
sample_posterior(state, nsamples = 1000L, thin = 10L, burnin = NULL,
  multiflip = TRUE, mcmc_args = list(), equilibrate_args = list(),
  history = FALSE, verbose = FALSE)
```

## Arguments

- state:

  A flat block state from
  [`block_state`](https://diogro.github.io/rgraphtool/reference/block_state.md),
  [`minimize_blockmodel_dl`](https://diogro.github.io/rgraphtool/reference/minimize_blockmodel_dl.md),
  [`weighted_block_state`](https://diogro.github.io/rgraphtool/reference/weighted_block_state.md),
  [`minimize_weighted_blockmodel_dl`](https://diogro.github.io/rgraphtool/reference/weighted_block_state.md),
  [`layered_block_state`](https://diogro.github.io/rgraphtool/reference/layered_block_state.md),
  [`minimize_layered_blockmodel_dl`](https://diogro.github.io/rgraphtool/reference/layered_block_state.md),
  [`pp_block_state`](https://diogro.github.io/rgraphtool/reference/pp_block_state.md),
  [`minimize_pp_blockmodel_dl`](https://diogro.github.io/rgraphtool/reference/pp_block_state.md),
  [`ranked_block_state`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md),
  or
  [`minimize_ranked_blockmodel_dl`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md),
  or a nested block state from
  [`nested_block_state`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md),
  [`minimize_nested_blockmodel_dl`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md),
  [`nested_weighted_block_state`](https://diogro.github.io/rgraphtool/reference/weighted_block_state.md),
  [`minimize_nested_weighted_blockmodel_dl`](https://diogro.github.io/rgraphtool/reference/weighted_block_state.md),
  [`nested_layered_block_state`](https://diogro.github.io/rgraphtool/reference/layered_block_state.md),
  [`minimize_nested_layered_blockmodel_dl`](https://diogro.github.io/rgraphtool/reference/layered_block_state.md),
  [`nested_ranked_block_state`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md),
  or
  [`minimize_nested_ranked_blockmodel_dl`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md).

- beta:

  Inverse temperature (\\\beta \ge 0\\, or `Inf` for greedy
  zero-temperature moves). Posterior sampling uses `beta = 1`.

- niter:

  Number of sweeps per call (for `mcmc_sweep`), or number of exponential
  cooling steps (for `mcmc_anneal`).

- c:

  Ergodic proposal smoothing parameter (\\c \ge 0\\; `Inf` gives uniform
  proposals). Defaults to `0.5` for flat states and `1` (scaled by
  \\2^l\\ at level \\l\\) for nested states. Nested states also accept a
  numeric vector of per-level values.

- d:

  Probability in \\\[0, 1)\\ of proposing a new (empty) group during a
  single-node move.

- multiflip:

  If `TRUE` (the default), use merge-split MCMC combining single-node
  moves, group splits, group merges, and merge-split moves; if `FALSE`,
  use single-node Metropolis-Hastings sweeps only.

- B_min, B_max:

  Bounds on the number of occupied groups (applied to level 0 for nested
  states).

- psingle, psplit, pmerge, pmergesplit:

  Relative move-type weights for merge-split MCMC. `psingle = NULL`
  defaults to the number of nodes at the current level.

- gibbs_sweeps:

  Number of Gibbs sweeps used to refine split and merge-split proposals
  (at least 1).

- ls_shuffle:

  For nested states, whether to randomize the order in which hierarchy
  levels are swept.

- wait, nbreaks, max_niter, force_niter, epsilon:

  Equilibration controls matching graph-tool's `mcmc_equilibrate()`:
  waits for `nbreaks` consecutive intervals of `wait` iterations without
  a record-breaking change in description length, up to `max_niter`
  iterations (or runs exactly `force_niter` iterations when specified).

- mcmc_args:

  Named list of arguments forwarded to `mcmc_sweep()`.

- history:

  Whether to return a data frame recording `nattempts`, `nmoves`,
  `entropy`, and `beta` at each iteration.

- callback:

  Optional R function called as `callback(state)` after each iteration.
  When `NULL`, the entire equilibration loop runs in C++.

- collect:

  Whether to collect sampled partitions (every `thin` iterations) during
  equilibration.

- thin:

  Thinning interval (number of equilibration steps between collected
  samples).

- verbose:

  Whether to print iteration progress.

- beta_range:

  Numeric vector of length 2 giving the starting and ending inverse
  temperatures for exponential cooling.

- mcmc_equilibrate_args:

  Named list of equilibration arguments used at each annealing step.

- nsamples:

  Number of thinned posterior partitions to collect.

- burnin:

  Number of burn-in sweeps before sampling. When `NULL` (the default),
  runs `mcmc_equilibrate()` with `equilibrate_args` first.

- equilibrate_args:

  Named list of equilibration controls used when `burnin = NULL`.

## Details

Single-node moves use Peixoto's neighbor-based ergodic proposal with
exact forward and reverse Hastings probabilities. Merge-split sweeps
combine single-node moves with group splits, group merges, and
merge-split moves refined by `gibbs_sweeps` conditional Gibbs sweeps
using an auxiliary-variable Hastings ratio that preserves exact detailed
balance for any `gibbs_sweeps >= 1` and any finite \\\beta \ge 0\\.

For nested states at level \\l \< L - 1\\, new groups created by
single-node moves or splits inherit the parent label at level \\l + 1\\;
correspondingly, emptying a singleton group or merging two groups at
level \\l\\ is restricted to groups sharing the same parent, while moves
at level \\l + 1\\ reassign groups across parents. This maintains
ergodicity and exact detailed balance over unlabelled hierarchies.

Unlike Python's mutating methods, these R functions treat states as
immutable and return the updated state in the `$state` element of the
result. All random draws use R's RNG, so
[`set.seed()`](https://rdrr.io/r/base/Random.html) reproduces
trajectories. As with any MCMC sampler on multimodal posteriors,
record-breaking equilibration is a heuristic and does not guarantee
global mixing across well-separated modes.

## Value

`mcmc_sweep()`, `mcmc_equilibrate()`, and `mcmc_anneal()` return a list
with the updated `state`, final `entropy`, cumulative `dS`, `nattempts`,
`nmoves`, `niter`, and optionally `history` and `samples`.
`sample_posterior()` returns a list with final `state`, `samples` (list
of 1-based partitions or nested hierarchies), `nested`, `equilibration`,
and `sampling`.

## References

Peixoto, T. P. (2014). Efficient Monte Carlo and greedy heuristic for
the inference of stochastic block models. *Physical Review E* 89,
012804. doi:10.1103/PhysRevE.89.012804

Peixoto, T. P. (2020). Merge-split Markov chain Monte Carlo for
community detection. *Physical Review E* 102, 012305.
doi:10.1103/PhysRevE.102.012305

## See also

[`partition_mode_state`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md),
[`mode_cluster_state`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md),
[`plot_marginals`](https://diogro.github.io/rgraphtool/reference/plot_marginals.md),
`demo("posterior_sampling")`

## Examples

``` r
g <- graph_create(12, directed = FALSE)
graph_add_edges(g, rbind(t(combn(1:6, 2)), t(combn(7:12, 2)), c(1, 7)))
set.seed(1)
fit <- minimize_blockmodel_dl(g)
post <- sample_posterior(fit, nsamples = 50, thin = 2, burnin = 20)
pm <- partition_mode_state(post)
pm
#> <partition mode: M = 50, B = 3, uncertainty = 0.1200, DL = 257.5850 nats>
```
