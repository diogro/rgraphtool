# Force-directed SFDP layout and breadth-first search traversal

`sfdp_layout()` computes a two-dimensional spring-electrical
force-directed layout in native C++ following Hu (2005) and graph-tool's
`sfdp_layout()`, using exact pairwise repulsion for small graphs (\\N
\le 128\\) and a Barnes-Hut quadtree for larger graphs.

`graph_bfs()` computes a breadth-first search traversal and
shortest-path spanning tree from `source` up to depth `max_dist`.

## Usage

``` r
sfdp_layout(graph, pos = NULL, C = 0.2, K = 1.0, p = 2.0, theta = 0.6,
  gamma = 1.0, mu = 0.0, mu_p = 1.0, init_step = NULL, cooling_step = 0.95,
  max_level = 10L, tol = 1e-4, max_iter = 200L, simple = FALSE,
  groups = NULL, pin = NULL)

graph_bfs(graph, source, max_dist = Inf, directed = NULL)
```

## Arguments

- graph:

  A graph created by
  [`graph_create`](https://diogro.github.io/rgraphtool/reference/graph.md).

- pos:

  Optional \\N \times 2\\ numeric matrix of initial vertex coordinates.
  Useful for warm-starting or incremental layout updates
  (`max_iter = 1L`).

- C:

  Relative repulsive force strength (positive scalar).

- K:

  Optimal spring length (positive scalar).

- p:

  Repulsive force exponent (positive scalar; default `2.0`).

- theta:

  Barnes-Hut quadtree opening angle parameter (used when \\N \> 128\\
  and `max_level > 0`).

- gamma:

  Inter-group force coefficient (`> 0` attracts distinct groups toward
  each other, `< 0` repels them).

- mu:

  Intra-group attractive force coefficient (`>= 0`).

- mu_p:

  Intra-group attractive scaling exponent.

- init_step:

  Optional positive initial displacement step size. Defaults to
  `max(K, 0.1 * K * sqrt(N))`.

- cooling_step:

  Multiplicative step-size cooling factor in `(0, 1]`.

- max_level:

  Maximum Barnes-Hut quadtree depth (integer `>= 0`).

- tol:

  Stopping tolerance on mean vertex displacement relative to `K`.

- max_iter:

  Maximum number of layout iterations (`0L` means no cap).

- simple:

  Logical; if `TRUE`, use monotone cooling
  (`step <- step * cooling_step`); if `FALSE`, use Hu (2005) adaptive
  step sizing.

- groups:

  Optional vector of length \\N\\ assigning vertices to groups for
  inter- and intra-group layout forces.

- pin:

  Optional logical vector of length \\N\\; vertices with `TRUE` remain
  fixed at their initial `pos` coordinates.

- source:

  1-based integer index of the starting vertex for BFS.

- max_dist:

  Maximum BFS depth (`Inf` for unrestricted traversal, or a non-negative
  integer).

- directed:

  Logical or `NULL`. If `NULL`, defaults to `isTRUE(graph$directed)`.

## Details

Both functions are deterministic given R's random number seed
([`set.seed`](https://rdrr.io/r/base/Random.html) for `sfdp_layout()`,
and completely deterministic for `graph_bfs()`). Self-loops are ignored
during layout and BFS traversal; parallel edges are collapsed to a
single tree edge in `graph_bfs()`.

## Value

`sfdp_layout()` returns an \\N \times 2\\ numeric matrix with column
names `"x"` and `"y"`.

`graph_bfs()` returns an object of class `"rgraphtool_bfs"` with
elements `source`, `max_dist`, `directed`, `distance` (numeric vector of
length \\N\\, `Inf` for unreached vertices), `parent` (integer vector of
1-based predecessors, `NA_integer_` for the source and unreached
vertices), `order` (1-based discovery order), and `tree_edges` (logical
vector aligned with
[`graph_edges`](https://diogro.github.io/rgraphtool/reference/graph.md)).

## Examples

``` r
g <- graph_create(10L, directed = FALSE)
graph_add_edges(g, cbind(1:9, 2:10))
set.seed(1)
pos <- sfdp_layout(g, max_iter = 40L)
bfs <- graph_bfs(g, source = 1L, max_dist = 3L)
print(bfs)
#> rgraphtool BFS traversal
#>   source: 1 (undirected, max_dist = 3)
#>   reached: 4 / 10 vertices (max distance = 3, tree edges = 3)
```
