# Create, modify, and inspect a graph and its property maps

Create and modify a native graph-tool adjacency structure, inspect its
vertices and edges, attach numeric edge property maps, and transform
edge weights with exact change-of-variables log-Jacobians.

## Usage

``` r
graph_create(n = 0L, directed = TRUE)
graph_add_vertices(graph, n = 1L)
graph_add_edges(graph, edges, properties = NULL)
graph_num_vertices(graph)
graph_num_edges(graph)
graph_edges(graph)
graph_degree(graph, mode = c("total", "in", "out"))
graph_set_edge_property(graph, name, values)
graph_edge_property(graph, name)
graph_edge_properties(graph)
graph_set_vertex_property(graph, name, values)
graph_vertex_property(graph, name)
graph_vertex_properties(graph)
graph_set_graph_property(graph, name, value)
graph_graph_property(graph, name)
graph_graph_properties(graph)
# S3 method for class 'rgraphtool_graph'
print(x, max_edges = 20L, max_ascii = 12L,
  ascii = NULL, ...)
transform_weights(x, transform = c("identity", "log", "log1p", "sqrt",
  "sign", "abs", "rank", "zscore"), shift = 0)
```

## Arguments

- n:

  Non-negative number of vertices.

- directed:

  Whether the graph is directed.

- graph:

  A graph created by `graph_create()`.

- edges:

  Two-column numeric matrix of 1-based source and target vertex IDs.

- properties:

  Optional named list or data frame of edge vectors (each of length
  `nrow(edges)`) to attach as edge property maps aligned with the newly
  added edges.

- mode:

  Degree type. Undirected graphs support only total degree.

- name:

  Non-empty character string naming a graph, vertex, or edge property
  map.

- values:

  Vector, matrix, or list of property values aligned with
  `graph_edges(graph)` (for edge properties) or
  `1:graph_num_vertices(graph)` (for vertex properties).

- value:

  Single value (scalar, vector, or raw object) for a graph-level
  property.

- x:

  For [`print()`](https://rdrr.io/r/base/print.html), an
  `rgraphtool_graph` object; for `transform_weights()`, a finite numeric
  vector of edge weights to transform.

- max_edges:

  Maximum number of edges to display in the edge preview.

- max_ascii:

  Maximum vertex count for automatically displaying the ASCII adjacency
  matrix.

- ascii:

  Optional logical overriding whether to display the ASCII adjacency
  matrix.

- ...:

  Unused; included for S3 `print` compatibility.

- transform:

  Weight transformation to apply. Invertible continuous transformations
  (`"identity"`, `"log"`, `"log1p"`, `"sqrt"`, `"zscore"`) attach the
  exact log-Jacobian \\\sum_e \ln \|f'(x_e)\|\\ as attribute
  `"log_jacobian"` for change-of-variables model comparison.

- shift:

  Optional finite numeric constant added to `x` before applying
  `transform`.

## Details

Graphs are mutable and aliases share the same underlying graph. External
pointers cannot be restored from R serialization. Parallel edges and
self-loops are allowed; a self-loop contributes two to total degree.
Edge endpoints are returned in their original orientation even for
undirected graphs.

Edge properties passed through `graph_add_edges(properties = ...)` are
stored by internal edge insertion index and returned by
`graph_edge_property()` in the exact row order of `graph_edges(graph)`,
even when edges are inserted out of vertex order.

## Value

`graph_create()` returns a graph. Modification functions invisibly
return the graph. Counts return numbers, `graph_edges()` returns a
two-column integer matrix, `graph_degree()` and `graph_edge_property()`
return numeric vectors, `graph_edge_properties()` returns a character
vector of property names, and `transform_weights()` returns a numeric
vector with attributes `"transform"` and `"log_jacobian"`.

## Examples

``` r
g <- graph_create(3)
graph_add_edges(g, rbind(c(1, 2), c(2, 3)),
                properties = list(weight = c(1.5, 4.2)))
graph_degree(g, "out")
#> [1] 1 1 0
w <- graph_edge_property(g, "weight")
transform_weights(w, "log")
#> [1] 0.4054651 1.4350845
#> attr(,"transform")
#> [1] "log"
#> attr(,"log_jacobian")
#> [1] -1.84055
```
