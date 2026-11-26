# Shortest paths, connected components, clustering, and assortativity

Native C++23 graph topology and descriptive network algorithms:
unweighted (BFS), non-negative weighted (Dijkstra), and negative-weight
(Bellman–Ford) shortest distances and path reconstruction, exact and
pseudo-diameter, connected and strongly connected components (with
attractor detection), Hopcroft–Tarjan biconnected components
(articulation points and bridges), topological sorting and DAG testing,
local and global clustering coefficients (with jackknife standard
error), and scalar/categorical assortativity.

## Usage

``` r
shortest_distance(graph, source = NULL, target = NULL, weights = NULL,
                  max_dist = Inf)
shortest_path(graph, source, target, weights = NULL)
pseudo_diameter(graph, source = 1L, weights = NULL)
graph_diameter(graph, weights = NULL)
label_components(graph)
label_biconnected_components(graph)
topological_sort(graph)
is_dag(graph)
local_clustering(graph)
global_clustering(graph)
assortativity(graph, values = "total", categorical = FALSE, weights = NULL)
```

## Arguments

- graph:

  An `rgraphtool_graph` object.

- source:

  Optional 1-based source vertex ID(s). When `NULL` in
  `shortest_distance()`, computes all-pairs distances in parallel.

- target:

  Optional 1-based target vertex ID(s).

- weights:

  Optional edge property name or numeric vector of edge lengths /
  weights.

- max_dist:

  Maximum distance cutoff for `shortest_distance()`.

- values:

  For `assortativity()`, `"total"`, `"in"`, or `"out"` for degree
  assortativity, a vertex property name, or a length-\\(N\\) vector.

- categorical:

  Logical; if `TRUE`, computes discrete Newman assortativity rather than
  scalar Pearson correlation.

## Value

`shortest_distance()` returns a numeric vector (for a single `source`)
or matrix. `shortest_path()` returns `list(vertices, edges, distance)`.
`pseudo_diameter()` returns `list(diameter, endpoints)` and
`graph_diameter()` returns a scalar. `label_components()` returns
`list(membership, sizes, num_components, attractors)`.
`label_biconnected_components()` returns
`list(edge_component, num_components, articulation_points, bridges)`.
`topological_sort()` returns an integer vertex ordering and `is_dag()`
returns a logical scalar. `local_clustering()` returns a length-\\(N\\)
vector, `global_clustering()` returns a named length-2 vector
`c(clustering, std_err)`, and `assortativity()` returns a numeric
correlation coefficient.

## Examples

``` r
g <- graph_create(6L, directed = FALSE)
graph_add_edges(g, rbind(c(1, 2), c(2, 3), c(3, 1), c(3, 4), c(4, 5), c(5, 6), c(6, 4)))
shortest_path(g, 1L, 6L)
#> $vertices
#> [1] 1 3 4 6
#> 
#> $edges
#> [1] 3 4 7
#> 
#> $distance
#> [1] 3
#> 
label_biconnected_components(g)
#> $edge_component
#> [1] 3 3 3 2 1 1 1
#> 
#> $num_components
#> [1] 3
#> 
#> $articulation_points
#> [1] 3 4
#> 
#> $bridges
#> [1] 4
#> 
global_clustering(g)
#> clustering    std_err 
#>  0.6000000  0.1844551 
assortativity(g, "total")
#> [1] -0.1666667
```
