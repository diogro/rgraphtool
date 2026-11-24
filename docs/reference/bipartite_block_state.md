# Native bipartite network representation and Stochastic Block Model inference

Construct, validate, infer, and summarize flat, hierarchical (nested),
and weighted microcanonical Bipartite Stochastic Block Models (biSBM) in
native multi-threaded C++23.

## Usage

``` r
bipartite_graph(x, type = NULL, directed = FALSE, dim_labels = NULL,
  weight_mode = c("auto", "binary", "multigraph", "weight"),
  weight_name = "weight", covariates = NULL)
is_bipartite(graph, type = NULL)
as_incidence_matrix(graph, type = NULL, weight = NULL)
bipartite_block_state(graph, type = NULL, blocks = NULL,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  edge_prior = c("bipartite", "unipartite"), weights = NULL, types = NULL,
  params = NULL, microcanonical = FALSE, dim_labels = NULL)
bipartite_description_length(state)
minimize_bipartite_blockmodel_dl(graph, type = NULL,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  edge_prior = c("bipartite", "unipartite"), weights = NULL, types = NULL,
  params = NULL, microcanonical = FALSE, B_row_range = NULL, B_col_range = NULL,
  nruns = 1L, refine = TRUE, control = list(), init_blocks = NULL,
  dim_labels = NULL,
  weight_mode = c("auto", "binary", "multigraph", "weight"))
nested_bipartite_block_state(graph, type = NULL, levels = NULL,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  edge_prior = c("bipartite", "unipartite"), weights = NULL, types = NULL,
  params = NULL, microcanonical = FALSE, dim_labels = NULL)
nested_bipartite_description_length(state)
minimize_nested_bipartite_blockmodel_dl(graph, type = NULL,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  edge_prior = c("bipartite", "unipartite"), weights = NULL, types = NULL,
  params = NULL, microcanonical = FALSE, nruns = 1L, control = list(),
  dim_labels = NULL,
  weight_mode = c("auto", "binary", "multigraph", "weight"))
bipartite_membership(state, level = 1L)
bipartite_edge_counts(state, level = 1L)
bipartite_block_matrix(state, level = 1L,
  normalize = c("density", "counts", "mean_weight"))
compare_bipartite_models(graph, type = NULL, weights = NULL, types = NULL,
  params = NULL, microcanonical = FALSE,
  edge_prior = c("bipartite", "unipartite"), nested = TRUE, nruns = 3L,
  dim_labels = NULL)
# S3 method for class 'rgraphtool_bipartite_block_state'
print(x, ...)
# S3 method for class 'rgraphtool_nested_bipartite_block_state'
print(x, ...)
# S3 method for class 'rgraphtool_bipartite_model_comparison'
print(x, ...)
```

## Arguments

- x:

  For `bipartite_graph()`, a rectangular numeric/logical incidence
  matrix (or sparse `Matrix`), an `rgraphtool_graph`, or an `igraph`
  object with a `type` vertex attribute. For
  [`print()`](https://rdrr.io/r/base/print.html), a bipartite state or
  comparison object.

- graph:

  An `rgraphtool_graph`, rectangular incidence matrix, or `igraph`
  object.

- type:

  Optional length-\\N\\ vertex type vector (logical, 0/1, 1/2,
  character, or factor with at most two distinct values). If `NULL`,
  reads vertex property `"type"` from `graph`.

- directed:

  Logical; whether the constructed bipartite graph is directed (default
  `FALSE`).

- dim_labels:

  Optional length-2 character vector naming the row and column vertex
  modes (e.g., `c("plant", "pollinator")`).

- weight_mode:

  How nonzero entries of an incidence matrix are encoded: `"binary"`
  (presence/absence), `"multigraph"` (non-negative integer parallel
  edges), `"weight"` (single edges with numeric edge property
  `weight_name`), or `"auto"`.

- weight_name:

  Edge property name used when `weight_mode = "weight"`.

- covariates:

  Optional named list of matrices with `dim(x)` to attach as additional
  edge properties.

- weight:

  Optional edge property name to aggregate in `as_incidence_matrix()`.

- blocks:

  Optional length-\\N\\ positive integer block assignment respecting
  bipartite vertex types.

- levels:

  Optional list of per-level integer block vectors for
  `nested_bipartite_block_state()`.

- degree_corrected:

  Logical; whether to fit the degree-corrected (DC) or
  non-degree-corrected (NDC) microcanonical bipartite SBM.

- degree_prior:

  Degree sequence hyperprior (`"distributed"` or `"uniform"`).

- edge_prior:

  Block edge-count matrix prior: `"bipartite"` (default, exact Yen &
  Larremore 2020 \\\ln\binom{B_1 B_2 + E - 1}{E}\\ prior exploiting
  bipartite structural zeros) or `"unipartite"` (matches `graph-tool`
  3.9 `BlockState(..., pclabel=type)`).

- weights:

  Optional edge weight property name(s), vector, or matrix for weighted
  bipartite SBMs.

- types:

  Optional covariate distribution type(s) (`"real-exponential"`,
  `"real-normal"`, `"discrete-geometric"`, `"discrete-poisson"`,
  `"discrete-binomial"`).

- params:

  Optional list of conjugate hyperparameter lists.

- microcanonical:

  Logical; whether to use microcanonical covariate priors.

- B_row_range, B_col_range:

  Optional length-2 integer bounds `c(min, max)` on the number of row
  and column blocks.

- nruns:

  Number of independent restarts; keeps the lowest description length.

- refine:

  Logical; whether to perform zero-temperature refinement sweeps.

- control:

  Optional list of search controls (`r` shrink ratio, `n_cand` merge
  candidates per block).

- init_blocks:

  Optional warm-start block assignment vector.

- state:

  A bipartite or nested bipartite block state.

- level:

  Hierarchy level (1-based) for nested bipartite states.

- normalize:

  Summary returned by `bipartite_block_matrix()`: `"density"` (\\e\_{rs}
  / (n_r n_s)\\), `"counts"` (\\e\_{rs}\\), or `"mean_weight"`.

- nested:

  Logical; whether `compare_bipartite_models()` also fits nested
  bipartite models.

- ...:

  Unused; included for S3 `print` compatibility.

## Details

The native C++ engine implements the exact microcanonical bipartite
description length of Yen & Larremore (2020, *Phys. Rev. E* 102, 032309)
with type-constrained partitions, exact per-type partition priors, and
the bipartite edge-count prior \\\ln\binom{B_1 B_2 + E - 1}{E}\\ (or
\\\ln\binom{2 B_1 B_2 + E - 1}{E}\\ for directed bipartite networks).
Setting `edge_prior = "unipartite"` reproduces `graph-tool` 3.9's
`BlockState(..., pclabel=type)` to machine precision.

Flat bipartite states also integrate with
[`mcmc_sweep()`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md),
[`mcmc_equilibrate()`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md),
and
[`sample_posterior()`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md)
for exact finite-temperature Bayesian posterior sampling.

## Value

`bipartite_graph()` returns an `rgraphtool_graph` with vertex property
`"type"`. `is_bipartite()` returns a logical scalar.
`as_incidence_matrix()`, `bipartite_edge_counts()`, and
`bipartite_block_matrix()` return rectangular numeric matrices.
`bipartite_block_state()` and `minimize_bipartite_blockmodel_dl()`
return an `rgraphtool_bipartite_block_state`; their nested counterparts
return an `rgraphtool_nested_bipartite_block_state`.
`bipartite_membership()` returns a named length-2 list of integer block
memberships for the two modes.

## Examples

``` r
set.seed(1)
M <- matrix(0L, 12, 12)
M[1:6, 1:6] <- rbinom(36, 1, 0.85)
M[7:12, 7:12] <- rbinom(36, 1, 0.85)
M[1:6, 7:12] <- rbinom(36, 1, 0.05)
M[7:12, 1:6] <- rbinom(36, 1, 0.05)

fit <- minimize_bipartite_blockmodel_dl(M, dim_labels = c("plant", "pollinator"))
fit
#> <bipartite block state (degree-corrected, unweighted, edge_prior = bipartite)>
#>   plant: 12 vertices -> 1 blocks
#>   pollinator: 12 vertices -> 1 blocks
#>   edges: 68 | description length: 132.3837 nats
bipartite_membership(fit)
#> $plant
#>  [1] 1 1 1 1 1 1 1 1 1 1 1 1
#> 
#> $pollinator
#>  [1] 1 1 1 1 1 1 1 1 1 1 1 1
#> 
bipartite_block_matrix(fit, normalize = "density")
#>          pollinator
#> plant     pollinator_1
#>   plant_1    0.4722222
```
