# Partition label alignment, vertex marginals, consensus, and posterior modes

Label alignment, vertex marginal distributions, MAP consensus
partitions, posterior uncertainty, and multimodal clustering of
posterior partition ensembles using Peixoto's random label model (2021).
These are the R counterparts of graph-tool's `align_partition_labels()`,
`partition_overlap()`, `PartitionModeState`, and `ModeClusterState`.

## Usage

``` r
align_partition_labels(x, y)
partition_overlap(x, y, norm = TRUE)
align_nested_partition_labels(x, y)
nested_partition_overlap(x, y, norm = TRUE)
partition_mode_state(partitions, relabel = TRUE, nested = NULL,
  converge = TRUE)
mode_marginals(mode, probability = FALSE, level = 0L)
mode_consensus(mode)
mode_uncertainty(mode, mle = TRUE)
mode_entropy(mode, mle = TRUE)
mode_lprob(mode, partition, mle = TRUE, align = FALSE)
mode_cluster_state(partitions, nested = NULL, b = NULL, relabel = TRUE,
  niter = 20L, converge = TRUE)
cluster_lprob(cluster, partition, mode = NULL, mle = TRUE, align = FALSE)
mode_summary(cluster)
```

## Arguments

- x, y:

  Flat integer label vectors (1-based) or nested lists of label vectors.

- norm:

  Whether to return the normalized overlap in \\\[0, 1\]\\ or the raw
  number of matching nodes.

- partitions:

  A list of flat or nested partitions, or the result object returned by
  [`sample_posterior`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md).

- relabel:

  Whether to align partition labels via maximum-weight bipartite
  matching when adding partitions.

- nested:

  Logical indicating whether partitions are hierarchical, or `NULL` to
  infer from `partitions`.

- converge:

  Whether to iterate label realignment until the random label model
  description length converges.

- mode:

  A partition mode object returned by `partition_mode_state()` or stored
  inside `mode_cluster_state()$modes`. For `mode_entropy()`, a
  `rgraphtool_mode_cluster` object is also accepted. For
  `cluster_lprob()`, an optional 1-based mode index (or `NULL` for the
  full mode-mixture log-probability).

- partition:

  A flat integer label vector or nested list of label vectors whose
  log-probability is evaluated under the mode or mode cluster.

- align:

  If `TRUE`, align `partition` to the mode's consensus labels before
  evaluating the log-probability. Defaults to `FALSE`, matching
  graph-tool's `posterior_lprob()`.

- probability:

  If `TRUE`, return row-normalized marginal probabilities; otherwise
  return integer counts.

- level:

  Hierarchy level (0-based) for nested modes.

- mle:

  If `TRUE`, use maximum-likelihood marginal probabilities; if `FALSE`,
  use Dirichlet(1) posterior mean estimates.

- b:

  Optional initial cluster assignment (1-based integer vector of length
  `length(partitions)`).

- niter:

  Maximum number of split/sweep/merge optimization iterations. Set
  `niter = 0` with `b` to score a fixed clustering.

- cluster:

  A mode cluster object returned by `mode_cluster_state()`.

## Details

Pairwise label alignment solves the exact maximum-weight bipartite
matching problem on the contingency matrix using the Hungarian
algorithm. For ensembles, `partition_mode_state()` aligns all partitions
to a common labelling by minimizing the integrated description length of
the random label model, matching graph-tool 3.9's `PartitionModeState`
entropy, `posterior_entropy()`, and `posterior_cdev()` numerically.

When the posterior distribution is multimodal, `mode_cluster_state()`
clusters the sampled partitions into distinct modes by minimizing the
mixture random label model description length, sorts the discovered
modes by decreasing weight, and aligns secondary modes to the dominant
mode.

## Value

`partition_mode_state()` returns a `rgraphtool_partition_mode` list with
`M`, `B`, `entropy`, `posterior_entropy`, `posterior_cdev`, `marginals`,
`consensus`, and aligned `partitions`. `mode_cluster_state()` returns a
`rgraphtool_mode_cluster` list with total `entropy`, `clusters`, number
of modes `K`, effective number of modes `Be`, `weights`, and `modes`.
`mode_summary()` formats the modes as a data frame.

## References

Peixoto, T. P. (2021). Revealing consensus and dissensus between network
partitions. *Physical Review X* 11, 021003.
doi:10.1103/PhysRevX.11.021003

## See also

[`sample_posterior`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md),
[`plot_marginals`](https://diogro.github.io/rgraphtool/reference/plot_marginals.md),
`demo("posterior_sampling")`, `demo("hierarchical_sampling")`

## Examples

``` r
bs <- list(c(1, 1, 1, 2, 2, 2), c(2, 2, 2, 1, 1, 1), c(1, 1, 2, 2, 2, 3))
pm <- partition_mode_state(bs)
pm
#> <partition mode: M = 3, B = 3, uncertainty = 0.1111, DL = 17.8045 nats>
mode_marginals(pm, probability = TRUE)
#>           [,1]      [,2]      [,3]
#> [1,] 0.0000000 1.0000000 0.0000000
#> [2,] 0.0000000 1.0000000 0.0000000
#> [3,] 0.3333333 0.6666667 0.0000000
#> [4,] 1.0000000 0.0000000 0.0000000
#> [5,] 1.0000000 0.0000000 0.0000000
#> [6,] 0.6666667 0.0000000 0.3333333
mode_consensus(pm)
#> [1] 2 2 2 1 1 1
mode_uncertainty(pm)
#> [1] 0.1111111
```
