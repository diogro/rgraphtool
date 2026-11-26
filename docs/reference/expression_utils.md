# Gene-Expression Fixed-Effects Projection, Filtering, Shrinkage Covariance, and Module Eigengenes

Utilities for conditioning gene-expression matrices on fixed-effect
designs via exact QR/REML orthogonal projection, pre-filtering
informative genes, computing analytical shrinkage covariances, and
computing sign-aligned module eigengenes.

## Usage

``` r
project_fixed_effects(
  expr,
  design,
  col_data = NULL,
  transform = c("zscore", "nonparanormal", "none")
)

filter_expression_genes(
  expr,
  col_data = NULL,
  nuisance = NULL,
  signal = NULL,
  min_var = 0,
  max_padj = 1,
  top_n = NULL
)

module_eigengenes(expr, blocks, level = 0L, min_size = 2L)

shrinkage_covariance(expr, target = c("diagonal", "identity"))
```

## Arguments

- expr:

  Numeric matrix of genes (rows) by samples (columns).

- design:

  A formula evaluated on `col_data` or an `M x p` numeric design matrix.

- col_data:

  Optional `data.frame` of sample metadata.

- transform:

  One of `"zscore"`, `"nonparanormal"`, or `"none"`.

- nuisance:

  Optional nuisance formula or matrix.

- signal:

  Optional signal formula or matrix tested beyond `nuisance`.

- min_var:

  Minimum residual variance after removing `nuisance`.

- max_padj:

  Maximum Benjamini-Hochberg adjusted p-value for `signal`.

- top_n:

  Optional maximum number of genes to retain.

- blocks:

  Integer block vector or `rgraphtool` block state.

- level:

  Hierarchy level (`0L` = finest) for nested states.

- min_size:

  Minimum block size for eigengene calculation.

- target:

  Shrinkage target (`"diagonal"` or `"identity"`).

## Value

`project_fixed_effects()` returns an object of class
`"rgraphtool_expression_projection"`; `filter_expression_genes()`
returns a list with `expr`, `stats`, and `keep`; `module_eigengenes()`
returns a list with `eigengenes`, `pve`, and `block_sizes`;
`shrinkage_covariance()` returns a list with `cov`, `lambda`, and
`n_obs`.
