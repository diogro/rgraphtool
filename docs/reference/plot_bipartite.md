# Visualize a fitted bipartite stochastic block model

Produces base-R visualizations of fitted bipartite stochastic block
models: block-reordered incidence matrices, mesoscale block-interaction
diagrams, and two-column bipartite network layouts.

## Usage

``` r
plot_bipartite(state, type = c("matrix", "meso", "network"), level = 1L,
  weight = NULL, main = NULL, palette = NULL, block_border = "#222222", ...)
```

## Arguments

- state:

  An `rgraphtool_bipartite_block_state` or
  `rgraphtool_nested_bipartite_block_state`.

- type:

  Visualization type: `"matrix"` (block-reordered incidence or weight
  heatmap with block boundaries), `"meso"` (two-column mesoscale block
  interaction diagram), or `"network"` (two-column bipartite node-link
  diagram colored by inferred block).

- level:

  Hierarchy level (1-based) for nested bipartite states.

- weight:

  Optional edge property name to plot in `type = "matrix"`.

- main:

  Optional plot title.

- palette:

  Optional color vector for `type = "matrix"`.

- block_border:

  Border color for block separators in `type = "matrix"`.

- ...:

  Additional graphical arguments passed to the underlying base-R
  plotting function.

## Value

Invisibly returns the plotted summary data (the reordered matrix and
orders for `"matrix"`, the density/count matrices for `"meso"`, or the
2D vertex coordinates for `"network"`).

## Examples

``` r
set.seed(1)
M <- matrix(0L, 10, 10)
M[1:5, 1:5] <- rbinom(25, 1, 0.8)
M[6:10, 6:10] <- rbinom(25, 1, 0.8)
fit <- minimize_bipartite_blockmodel_dl(M, dim_labels = c("host", "parasite"))
plot_bipartite(fit, type = "matrix")

plot_bipartite(fit, type = "meso")
```
