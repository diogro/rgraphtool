# Compare weighted SBM covariate distributions and transformations

Declarative model selection across weighted stochastic block model
covariate distributions and weight transformations (such as `log-normal`
vs `real-exponential` and `real-normal`), automatically applying the
exact change-of-variables Jacobian adjustment so transformed models are
compared on the original weight scale.

## Usage

``` r
compare_weight_models(graph, weights, candidates = NULL, nested = TRUE,
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

- weights:

  Single edge weight covariate: an edge property name or numeric vector
  of length `graph_num_edges(graph)`.

- candidates:

  Candidate weight models to compare:

  - `NULL` (default): automatically selects discrete count models
    (`"discrete-geometric"`, `"discrete-poisson"`,
    `"discrete-binomial"`) for non-negative integer weights, continuous
    positive models (`"real-exponential"`, `"real-normal"`,
    `"log-normal"`) for strictly positive real weights, or
    microcanonical vs Bayesian `"real-normal"` for signed real weights.

  - A character vector of shorthand names from `"real-exponential"`,
    `"real-normal"`, `"log-normal"`, `"discrete-geometric"`,
    `"discrete-poisson"`, `"discrete-binomial"`.

  - A named list where each entry is a list with elements `type`, and
    optional `transform` (default `"identity"`), `shift` (default `0`),
    and `params` (default `"microcanonical"`).

- nested:

  Whether to fit nested (`TRUE`, default) or flat (`FALSE`) weighted
  SBMs.

- degree_corrected:

  Whether to use degree correction.

- degree_prior:

  Degree sequence hyperprior (`"distributed"` or `"uniform"`).

- method:

  Rank candidate models by `"map"` description length or full Bayesian
  `"posterior"` marginal evidence.

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

  Whether to run refinement sweeps during minimization.

## Details

For an invertible continuous weight transformation \\(y_e = f(w_e)\\),
the induced density on the original weights is \\(p_W(w) = p_Y(f(w))
\prod_e \|f'(w_e)\|\\), so the adjusted description length and negative
log-evidence are \$\$-\ln p_W(w) = -\ln p_Y(f(w)) - \sum_e \ln
\|f'(w_e)\|.\$\$ For `log-normal` (\\(f(w) = \ln w\\)), the Jacobian
adjustment added to the transformed description length is \\(+\sum_e \ln
w_e\\).

Note that discrete probability mass functions and continuous probability
density functions live on different reference measures; by default
(`candidates = NULL`), discrete weights are compared against discrete
models and continuous weights against continuous models.

## Value

A `rgraphtool_weight_model_comparison` list containing:

- `summary`: ranked data frame with `model`, `type`, `transform`, `B`,
  `raw_dl`, `jacobian_adj`, `description_length`, `neg_log_evidence`,
  `posterior_entropy`, `delta`, and `posterior_weight`.

- `best_model`: name of the winning candidate.

- `best_state`: fitted weighted state for the winning candidate.

- `states`: named list of fitted states ordered from best to worst.

- `evidences`: named list of `rgraphtool_model_evidence` objects (when
  `method = "posterior"`).

- `log_bayes_factors`: pairwise matrix of log Bayes factors.

## References

Peixoto, T. P. (2018). Nonparametric weighted stochastic block models.
*Physical Review E* 97, 012306.

## See also

[`weighted_block_state`](https://diogro.github.io/rgraphtool/reference/weighted_block_state.md),
[`transform_weights`](https://diogro.github.io/rgraphtool/reference/graph.md),
[`model_evidence`](https://diogro.github.io/rgraphtool/reference/compare_model_classes.md),
`demo("weighted_sbm")`

## Examples

``` r
set.seed(42)
g <- graph_create(12, directed = FALSE)
edges <- rbind(t(utils::combn(1:6, 2)), t(utils::combn(7:12, 2)), c(6, 7))
w <- c(exp(stats::rnorm(15, 2, 0.8)), exp(stats::rnorm(15, 0, 0.8)), 1.5)
graph_add_edges(g, edges, properties = list(w = w))
cmp <- compare_weight_models(g, weights = "w", nested = FALSE)
cmp
#> <weighted SBM model comparison (flat, DC, MAP)>
#>             model B description_length jacobian_adj  delta posterior_weight
#>        log-normal 2            120.890       32.052  0.000         9.99e-01
#>       real-normal 2            127.987        0.000  7.097         8.27e-04
#>  real-exponential 2            131.065        0.000 10.175         3.81e-05
#> Best weight model: log-normal
```
