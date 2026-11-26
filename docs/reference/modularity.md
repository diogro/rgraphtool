# Newman modularity and modularity maximization (use with extreme caution)

Computes generalized Newman–Girvan / Reichardt–Bornholdt modularity
\\(Q(\gamma)\\) and maximizes it via fast native C++23 multilevel
Louvain coarsening and fine-grained refinement.

## Usage

``` r
modularity(graph, blocks, gamma = 1.0, weights = NULL)
maximize_modularity(graph, gamma = 1.0, weights = NULL, nruns = 1L,
                    max_levels = 10L, warn = TRUE)
# S3 method for class 'rgraphtool_modularity_state'
print(x, ...)
```

## Arguments

- graph:

  An `rgraphtool_graph` object.

- blocks:

  Integer partition vector of length `graph_num_vertices(graph)` (or a
  fitted block/modularity state).

- gamma:

  Resolution parameter \\(\gamma\\) (default `1.0`).

- weights:

  Optional edge property name or numeric edge weight vector.

- nruns:

  Number of randomized multilevel Louvain restarts; keeps the highest
  modularity score.

- max_levels:

  Maximum number of Louvain coarsening passes.

- warn:

  Logical; if `TRUE` (the default), emits an explicit warning explaining
  why modularity maximization is statistically unsound.

- x:

  An `rgraphtool_modularity_state` object.

- ...:

  Unused; included for S3 `print` compatibility.

## DANGER: Modularity Maximization Considered Harmful

Using modularity maximization for community detection is almost always
**a terrible idea**. Following `graph-tool`'s upstream documentation,
`maximize_modularity()` is included *strictly for comparison and
pedagogical purposes* and emits an explicit warning by default:

- **Massive overfitting of random noise**: Modularity has no statistical
  regularization or prior penalty on the number of groups. On completely
  structureless Erdos–Renyi or configuration-model random graphs,
  `maximize_modularity()` routinely fabricates spurious partitions with
  high modularity (\\(Q \> 0.3\\)–\\(0.8\\)), whereas
  [`minimize_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/minimize_blockmodel_dl.md)
  rightly infers \\(B = 1\\).

- **Resolution limit (underfitting)**: Simultaneously, modularity cannot
  resolve genuine communities smaller than \\(\sim \sqrt{2E}\\) edges
  (Fortunato & Barthelemy 2007), merging distinct small modules in large
  networks.

- **Inferential alternative**: For principled Bayesian community
  detection—including specifically assortative communities—use
  [`minimize_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/minimize_blockmodel_dl.md),
  [`minimize_nested_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md),
  or
  [`minimize_pp_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/pp_block_state.md)
  instead.

See Tiago P. Peixoto (2023), *Descriptive vs. inferential community
detection in networks: pitfalls, myths and half-truths* (Cambridge
University Press) and <https://skewed.de/tiago/blog/modularity-harmful>.

## Value

`modularity()` returns the scalar modularity \\(Q(\gamma)\\).
`maximize_modularity()` returns an `rgraphtool_modularity_state` list
containing `blocks`, `modularity`, `B`, `gamma`, and `graph`.

## Examples

``` r
# Illustration of modularity overfitting a completely random Erdos-Renyi graph:
set.seed(42)
g_rand <- random_graph_erdos_renyi(80L, p = 0.06)

# Modularity maximization fabricates spurious groups (and warns!):
mod_fit <- suppressWarnings(maximize_modularity(g_rand))
mod_fit$B
#> [1] 7
mod_fit$modularity
#> [1] 0.4340741

# Principled MDL inference rightly finds B = 1 (no community structure):
sbm_fit <- minimize_blockmodel_dl(g_rand)
length(block_sizes(sbm_fit))
#> [1] 1
```
