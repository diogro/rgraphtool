# Analytical Ledoit-Wolf / Schafer-Strimmer Shrinkage Covariance

Computes an optimal linear shrinkage covariance matrix \\((1 - \lambda)
S + \lambda T\\) for high-dimensional data (\\(N \gg M\\)) with zero
external package dependencies, suitable for passing directly to
[`pseudo_normal_block_state`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md)`(cov = ..., n_eff = ...)`.

## Usage

``` r
shrinkage_covariance(x, target = c("diagonal", "identity"))
```

## Arguments

- x:

  Numeric `N x M` matrix of variables (rows) by observations (columns),
  with `M >= 3`.

- target:

  Shrinkage target: `"diagonal"` (unequal variable variances, zero
  covariances) or `"identity"` (mean-variance scaled identity).

## Value

A list with `cov` (`N x N` symmetric positive-definite matrix), `lambda`
(shrinkage intensity in `[0, 1]`), and `n_obs` (`M`).

## See also

[`pseudo_normal_block_state`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md)

## Examples

``` r
set.seed(1)
x <- matrix(rnorm(20 * 10), nrow = 20, ncol = 10)
sh <- shrinkage_covariance(x, target = "diagonal")
sh$lambda
#> [1] 0.7326759
dim(sh$cov)
#> [1] 20 20
```
