# Read, write, and convert graphs between graph-tool (.gt) and igraph

Native C++ reader and writer for the `graph-tool` binary `.gt` file
format (version 1, both little-endian and big-endian, with optional
gzip/bzip2/xz/zstd compression) and lossless bidirectional conversion
with igraph graph objects.

## Usage

``` r
read_gt(file)
write_gt(graph, file, comment = NULL,
         compression = c("auto", "none", "gzip", "bzip2", "xz", "zstd"))
graph_load(file, fmt = c("auto", "gt", "graphml", "gml", "edgelist"), ...)
graph_save(graph, file, fmt = c("auto", "gt", "graphml", "gml", "edgelist"), ...)
as_igraph(graph, include_properties = TRUE)
from_igraph(graph, include_properties = TRUE)
```

## Arguments

- file:

  Path to the input or output graph file.

- graph:

  An `rgraphtool_graph` object (for `write_gt()`, `graph_save()`, and
  `as_igraph()`) or an `igraph` object (for `from_igraph()`).

- comment:

  Optional single character string stored in the `.gt` header comment
  field. When `NULL`, a descriptive default header is generated.

- compression:

  Compression applied when writing `.gt` files. `"auto"` infers
  compression from the file extension (`.gz`, `.bz2`, `.xz`, `.zst`);
  `read_gt()` automatically detects compression from magic bytes.

- fmt:

  File format for `graph_load()` and `graph_save()`. `"auto"` infers the
  format from the file extension. `"gt"` uses the native C++ binary
  reader/writer without external dependencies; `"graphml"`, `"gml"`, and
  `"edgelist"` delegate to igraph.

- include_properties:

  Logical; whether to transfer graph, vertex, and edge property maps
  during igraph conversion.

- ...:

  Additional arguments forwarded to `write_gt()` (for `fmt = "gt"`) or
  to
  [`igraph::read_graph()`](https://r.igraph.org/reference/read_graph.html)
  /
  [`igraph::write_graph()`](https://r.igraph.org/reference/write_graph.html).

## Details

`read_gt()` and `write_gt()` implement the official `graph-tool` binary
`.gt` specification directly in C++ with no Python dependency. They
preserve directedness, parallel edges, self-loops, isolated vertices,
and all 15 `graph-tool` internal property map types (`bool`, `int16_t`,
`int32_t`, `int64_t`, `double`, `long double`, `string`, their
`vector<T>` counterparts, and raw pickled `python::object` payloads)
across graph, vertex, and edge property maps. Numeric edge properties
loaded from `.gt` or imported via `from_igraph()` are stored directly in
C++ so weighted and layered SBM inference functions and C++ extensions
can access them immediately.

`as_igraph()` and `from_igraph()` require the optional igraph package
and preserve vertex counts, directedness, parallel edges, self-loops,
and graph, vertex, and edge attributes.

## Value

`read_gt()`, `graph_load()`, and `from_igraph()` return an
`rgraphtool_graph` object. `as_igraph()` returns an `igraph` object.
`write_gt()` and `graph_save()` invisibly return `file`.

## Examples

``` r
g <- graph_create(4L, directed = TRUE)
graph_add_edges(g, rbind(c(1, 2), c(2, 3), c(3, 1), c(2, 4)),
                properties = list(weight = c(1.2, 0.8, 2.5, 1.0)))
graph_set_vertex_property(g, "name", c("a", "b", "c", "d"))
graph_set_graph_property(g, "title", "example")

f <- tempfile(fileext = ".gt.gz")
write_gt(g, f)
g2 <- read_gt(f)
graph_edge_property(g2, "weight")
#> [1] 1.2 0.8 1.0 2.5
graph_vertex_property(g2, "name")
#> [1] "a" "b" "c" "d"
unlink(f)

if (requireNamespace("igraph", quietly = TRUE)) {
  ig <- as_igraph(g)
  g3 <- from_igraph(ig)
  graph_num_edges(g3)
}
#> [1] 4
```
