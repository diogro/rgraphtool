# Vertex and edge centrality measures

Multi-threaded native C++23 implementations of classic network
centrality measures matching `graph-tool`'s `centrality` module:
PageRank, Brandes vertex and edge betweenness, closeness (standard and
harmonic), leading eigenvector centrality, Katz centrality, and
Kleinberg HITS hubs and authorities.

## Usage

``` r
pagerank(graph, damping = 0.85, pers = NULL, weights = NULL,
         tol = 1e-6, max_iter = 1000L)
betweenness(graph, weights = NULL, normalized = TRUE)
closeness(graph, weights = NULL, normalized = TRUE, harmonic = FALSE)
eigenvector_centrality(graph, weights = NULL, tol = 1e-6, max_iter = 1000L)
katz_centrality(graph, alpha = 0.01, beta = NULL, weights = NULL,
                normalized = TRUE, tol = 1e-6, max_iter = 1000L)
hits_centrality(graph, weights = NULL, tol = 1e-6, max_iter = 1000L)
```

## Arguments

- graph:

  An `rgraphtool_graph` object.

- damping:

  PageRank damping factor in \\(\[0, 1\]\\) (default `0.85`).

- pers:

  Optional non-negative personalization vector of length
  `graph_num_vertices(graph)`.

- weights:

  Optional edge property name or numeric vector of length
  `graph_num_edges(graph)`. For `betweenness()` and `closeness()`,
  weights are interpreted as path lengths (distances); for spectral/walk
  centralities, weights are interpreted as connection strengths.

- tol:

  Convergence tolerance for iterative power methods.

- max_iter:

  Maximum number of power iterations.

- normalized:

  Logical; whether to normalize centrality values.

- harmonic:

  Logical; if `TRUE`, computes harmonic closeness \\(\sum\_{v \neq u} 1
  / d(u, v)\\) (divided by \\(N - 1\\) when `normalized = TRUE`).

- alpha:

  Katz attenuation factor (should be smaller than \\(1 /
  \lambda\_{\max}\\)).

- beta:

  Optional exogenous vertex weight vector of length \\(N\\) for
  `katz_centrality()`.

## Value

`pagerank()`, `closeness()`, and `katz_centrality()` return numeric
vectors of length \\(N\\). `betweenness()` returns a list with elements
`vertex` (length \\(N\\)) and `edge` (length \\(E\\)).
`eigenvector_centrality()` returns `list(eigenvalue, vector)`, and
`hits_centrality()` returns `list(eigenvalue, authority, hub)`.

## Examples

``` r
g <- graph_create(6L, directed = FALSE)
graph_add_edges(g, rbind(c(1, 2), c(2, 3), c(3, 1), c(3, 4), c(4, 5), c(5, 6), c(6, 4)))
pagerank(g)
#> [1] 0.1459854 0.1459854 0.2080292 0.2080292 0.1459854 0.1459854
#> attr(,"iterations")
#> [1] 8
betweenness(g)
#> $vertex
#> [1] 0.0 0.0 0.6 0.6 0.0 0.0
#> 
#> $edge
#> [1] 0.06666667 0.26666667 0.26666667 0.60000000 0.26666667 0.06666667 0.26666667
#> 
closeness(g)
#> [1] 0.5000000 0.5000000 0.7142857 0.7142857 0.5000000 0.5000000
```
