# Latent triadic closure stochastic block model

Native C++ implementation of Peixoto's (2022) joint stochastic block
model with \\L \ge 1\\ generations of zero-inflated triadic closures,
matching graph-tool's `LatentClosureBlockState`. Disentangles observed
edges into homophily/community substrate edges versus triadic-closure
edges and estimates per-vertex open/closed triad counts and
transitivity.

## Usage

``` r
latent_closure_block_state(graph, L = 1L, nested = TRUE,
  blocks = NULL, gen = NULL, w0 = NULL, loops = NULL, ex = NULL,
  degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"), q_max = Inf)

closure_description_length(state, breakdown = FALSE)

infer_triadic_closure(state, nsamples = 200L, burnin = 100L, thin = 1L,
  beta = 1, pedges = 0.5, c = 1, d = 0.01,
  sample_loops = TRUE, init_sbm = TRUE,
  collect_partitions = FALSE)
```

## Arguments

- graph:

  An undirected simple graph created by
  [`graph_create`](https://diogro.github.io/rgraphtool/reference/graph.md)
  (no self-loops or parallel edges).

- L:

  Maximum number of triadic closure generations (\\L \ge 1\\).

- nested:

  Whether to use a nested (`TRUE`) or flat (`FALSE`) SBM substrate
  prior.

- blocks:

  Optional initial substrate partition or hierarchy.

- gen:

  Optional integer vector of length `graph_num_edges(graph)` in `0..L`
  giving each edge's initial generation (`0` = substrate, `1..L` =
  triadic closure generation). Defaults to all `0L`.

- w0:

  Optional initial substrate edge multiplicities for generation-0 edges.

- loops:

  Optional length-\\N\\ integer vector of initial substrate self-loop
  multiplicities.

- ex:

  Optional list of length `graph_num_edges(graph)` giving the 1-based
  closing ego vertex IDs for each closure edge (empty integer vectors
  for substrate edges).

- degree_corrected:

  Whether the substrate SBM is degree-corrected.

- degree_prior:

  Degree sequence hyperprior (`"distributed"` or `"uniform"`).

- q_max:

  Maximum degree cached in the partition-count table.

- state:

  A `rgraphtool_latent_closure_block_state` object.

- breakdown:

  If `TRUE`, return a named list with `total`, `sbm`, `closure`,
  per-generation `layer_entropies`, and per-vertex open/closed triad
  counts `M` and `E`.

- nsamples, burnin, thin:

  Posterior MCMC sample count, burn-in sweeps, and thinning interval.

- beta:

  Inverse MCMC temperature.

- pedges:

  Probability of performing an edge-layer/ego sweep versus a substrate
  SBM partition sweep.

- c, d:

  Substrate SBM proposal parameters.

- sample_loops:

  Whether to sample latent substrate self-loop multiplicities during
  edge sweeps.

- init_sbm:

  If `TRUE` and `state` starts from a single-group partition, initialize
  the substrate SBM partition before MCMC.

- collect_partitions:

  If `TRUE`, return sampled substrate SBM partitions in `partitions`.

## Value

`latent_closure_block_state()` returns a
`rgraphtool_latent_closure_block_state`. `closure_description_length()`
returns the exact joint description length (or a detailed breakdown
list). `infer_triadic_closure()` returns a
`rgraphtool_closure_inference` containing per-edge posterior layer
probabilities (`edges`: `substrate_prob`, `closure_prob`, and
`gen_0`..`gen_L`), per-generation vertex statistics (`vertex_stats`:
average closed triads `t`, open triads `m`, and transitivity fraction
`c`), the final `state`, and entropy trace.

## References

Peixoto, T. P. (2022). Disentangling homophily, community structure, and
triadic closure in networks. *Physical Review X* 12, 011004.
doi:10.1103/PhysRevX.12.011004

## See also

[`measured_block_state`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md),
[`predict_edges`](https://diogro.github.io/rgraphtool/reference/predict_edges.md),
`demo("triadic_closure")`

## Examples

``` r
g <- graph_create(6, directed = FALSE)
graph_add_edges(g, rbind(cbind(1:5, 2:6), c(6L, 1L), c(1L, 3L), c(3L, 5L)))
st <- latent_closure_block_state(g, L = 1L, nested = FALSE)
set.seed(1)
res <- infer_triadic_closure(st, nsamples = 40L, burnin = 30L, init_sbm = FALSE)
res$edges
#>   from to substrate_prob closure_prob gen_0 gen_1
#> 1    1  2          0.900        0.100 0.900 0.100
#> 2    1  3          0.900        0.100 0.900 0.100
#> 3    2  3          1.000        0.000 1.000 0.000
#> 4    3  4          1.000        0.000 1.000 0.000
#> 5    3  5          0.875        0.125 0.875 0.125
#> 6    4  5          1.000        0.000 1.000 0.000
#> 7    5  6          1.000        0.000 1.000 0.000
#> 8    1  6          1.000        0.000 1.000 0.000
```
