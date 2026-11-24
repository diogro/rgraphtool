# Uncertain network reconstruction and latent Poisson multigraphs

Native C++ implementations of Bayesian network reconstruction from noisy
or uncertain measurements (Peixoto 2018) and erased Poisson multigraph
inference (Peixoto 2020), matching graph-tool's `MeasuredBlockState`,
`MixedMeasuredBlockState`, `UncertainBlockState`,
`LatentMultigraphBlockState`, and `latent_multigraph()`.

## Usage

``` r
measured_block_state(graph, n = NULL, x = NULL,
  n_default = 1L, x_default = 0L,
  fn_params = list(alpha = 1, beta = 1),
  fp_params = list(mu = 1, nu = 1),
  lp = NA_real_, lq = NA_real_,
  nested = TRUE, blocks = NULL,
  latent_edges = NULL, latent_weights = NULL,
  degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"),
  self_loops = FALSE, q_max = Inf)

mixed_measured_block_state(graph, n = NULL, x = NULL,
  n_default = 1L, x_default = 0L,
  fn_params = list(alpha = 1, beta = 1),
  fp_params = list(mu = 1, nu = 1),
  nested = TRUE, blocks = NULL,
  latent_edges = NULL, latent_weights = NULL,
  degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"),
  self_loops = FALSE, q_max = Inf)

uncertain_block_state(graph, q = NULL, q_default = 0,
  nested = TRUE, blocks = NULL,
  latent_edges = NULL, latent_weights = NULL,
  degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"),
  self_loops = FALSE, q_max = Inf)

latent_multigraph_block_state(graph, nested = TRUE, blocks = NULL,
  latent_weights = NULL, degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"),
  self_loops = TRUE, q_max = Inf)

uncertain_description_length(state, sbm = TRUE, density = FALSE,
  aE = 1, breakdown = FALSE)

measured_error_posterior(state)

edge_probabilities(state, edges, log = TRUE, epsilon = 1e-8)

reconstruct_network(state, nsamples = 200L, burnin = 100L, thin = 1L,
  beta = 1, pedges = 0.5, ph = 0.1, hstep = 1,
  c = 1, d = 0.01, init_sbm = TRUE,
  collect_partitions = FALSE, threshold = 0.5)

latent_multigraph(graph, epsilon = 1e-8, max_niter = 0L)
```

## Arguments

- graph:

  A graph created by
  [`graph_create`](https://diogro.github.io/rgraphtool/reference/graph.md).

- n, x:

  Total measurement trials \\n\_{ij}\\ and positive observations
  \\x\_{ij}\\ for each observed edge of `graph`, supplied as a numeric
  vector of length `graph_num_edges(graph)`, an edge property name, or
  `NULL` (defaulting to 1).

- n_default, x_default:

  Default measurement counts for unobserved dyads not listed in `graph`.

- fn_params:

  Named list `list(alpha, beta)` giving Beta prior parameters for the
  false-negative error rate \\p\\.

- fp_params:

  Named list `list(mu, nu)` giving Beta prior parameters for the
  false-positive error rate \\q\\.

- lp, lq:

  Optional fixed log error rates \\\ln p\\ and \\\ln q\\ (\\\le 0\\).
  Leave as `NA_real_` to integrate out unknown error rates under their
  Beta priors.

- q:

  Extraneous edge probability estimates \\q\_{ij} \in \[0, 1\]\\ for
  observed edges of `graph`, given as a numeric vector, edge property
  name, or `NULL`.

- q_default:

  Default probability in \\\[0, 1\]\\ for unobserved dyads.

- nested:

  Whether to place a nested hierarchical SBM prior (`TRUE`) or a flat
  SBM prior (`FALSE`) on the latent multigraph.

- blocks:

  Optional initial partition: `NULL`, a length-\\N\\ integer vector, a
  hierarchy list, or a fitted block state.

- latent_edges, latent_weights:

  Optional initial latent multigraph edgelist (2-column matrix) and
  integer multiplicities. Defaults to the edges of `graph` with unit
  weight.

- degree_corrected:

  Whether the SBM prior is degree-corrected.

- degree_prior:

  Degree sequence hyperprior (`"distributed"` or `"uniform"`).

- self_loops:

  Whether self-loops are valid dyads in the latent graph.

- q_max:

  Maximum degree cached in the partition-count table.

- state:

  A `rgraphtool_uncertain_block_state` object.

- sbm:

  Include the SBM prior description length (`TRUE`) or return only the
  measurement negative log-likelihood (`FALSE`).

- density:

  Include an optional Poisson edge-count density prior with mean `aE`.

- aE:

  Expected total number of edges when `density = TRUE`.

- breakdown:

  If `TRUE`, return a named list separating total, measurement, and SBM
  description lengths plus Beta posterior parameters.

- edges:

  Two-column matrix or data.frame of 1-based vertex pairs to query.

- log:

  Return log-probabilities (`TRUE`) or probabilities (`FALSE`).

- epsilon:

  Convergence tolerance for latent multiplicity summation or EM
  iteration.

- nsamples, burnin, thin:

  Number of posterior samples collected, burn-in sweeps discarded, and
  sweeps between consecutive samples.

- beta:

  Inverse MCMC temperature (default 1 for posterior sampling).

- pedges:

  Probability of performing an edge-multiplicity sweep versus an SBM
  partition sweep.

- ph, hstep:

  Hyperparameter move probability and proposal step size for
  `mixed_measured_block_state()`.

- c, d:

  SBM proposal mixing and empty-group creation parameters.

- init_sbm:

  If `TRUE` and `state` starts from a single-group partition, initialize
  the SBM hierarchy on the initial latent graph before MCMC.

- collect_partitions:

  If `TRUE`, return sampled SBM partitions in `partitions`.

- threshold:

  Marginal edge-probability threshold used to construct the returned
  `consensus_graph`.

- max_niter:

  Maximum EM iterations for `latent_multigraph()` (`0L` means run until
  `epsilon` convergence).

## Details

All four reconstruction models combine an exact dynamic-edge flat or
nested SBM prior with a measurement likelihood over observed dyads:

- `measured_block_state()`: repeated binary observations \\(n\_{ij},
  x\_{ij})\\ with unknown Beta-prior or fixed homogeneous false-negative
  (\\p\\) and false-positive (\\q\\) rates.

- `mixed_measured_block_state()`: heterogeneous per-dyad Beta-Binomial
  error rates with joint hyperparameter MCMC over \\(\alpha, \beta, \mu,
  \nu)\\.

- `uncertain_block_state()`: extraneous per-edge probability estimates
  \\q\_{ij}\\ and default non-edge probability
  \\q\_{\mathrm{default}}\\.

- `latent_multigraph_block_state()`: Bayesian posterior over latent edge
  multiplicities and self-loops conditioned on an erased simple graph.

`latent_multigraph()` implements the fast \\O(V + E)\\
Expectation-Maximization algorithm of Peixoto (2020) for
maximum-likelihood expected edge and self-loop multiplicities.

## Value

Constructors return a `rgraphtool_uncertain_block_state`.
`uncertain_description_length()` returns a numeric scalar (or breakdown
list). `measured_error_posterior()` returns Beta posterior parameters,
means, and standard deviations for \\p\\ and \\q\\.
`edge_probabilities()` returns a numeric vector of conditional (log-)
probabilities. `reconstruct_network()` returns a
`rgraphtool_reconstruction` containing posterior marginal edge
probabilities (`edges`), a thresholded `consensus_graph`, the final
`state`, posterior error-rate traces, and entropy history.
`latent_multigraph()` returns a `rgraphtool_latent_multigraph` with
expected edge and self-loop weights (`edges`), vertex propensities
`theta_out`/`theta_in`, and convergence diagnostics.

## References

Peixoto, T. P. (2018). Reconstructing networks with unknown and
heterogeneous errors. *Physical Review X* 8, 041011.
doi:10.1103/PhysRevX.8.041011

Peixoto, T. P. (2020). Latent Poisson models for networks with
heterogeneous density. *Physical Review E* 102, 012309.
doi:10.1103/PhysRevE.102.012309

## See also

[`predict_edges`](https://diogro.github.io/rgraphtool/reference/predict_edges.md),
[`evaluate_edge_prediction`](https://diogro.github.io/rgraphtool/reference/predict_edges.md),
[`latent_closure_block_state`](https://diogro.github.io/rgraphtool/reference/latent_closure_block_state.md),
`demo("uncertain_reconstruction")`

## Examples

``` r
g <- graph_create(8, directed = FALSE)
e <- rbind(t(combn(1:4, 2)), t(combn(5:8, 2)), c(4, 5))
graph_add_edges(g, e, properties = list(n = rep(3, nrow(e)), x = rep(3, nrow(e))))
st <- measured_block_state(g, n = "n", x = "x", n_default = 3L, x_default = 0L)
uncertain_description_length(st)
#> [1] 32.68894
em <- latent_multigraph(g)
head(em$edges)
#>   from to   weight is_self_loop
#> 1    1  2 1.253635        FALSE
#> 2    1  3 1.253635        FALSE
#> 3    1  4 1.456230        FALSE
#> 4    2  3 1.253635        FALSE
#> 5    2  4 1.456230        FALSE
#> 6    3  4 1.456230        FALSE
```
