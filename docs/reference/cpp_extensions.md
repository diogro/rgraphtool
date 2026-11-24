# Native C++ extension interface, k-core decomposition, and property map access

Exposes rgraphtool's public C++ extension header (`<rgraphtool.h>`) for
writing high-performance custom network algorithms in C++ via
[`sourceCpp`](https://rdrr.io/pkg/Rcpp/man/sourceCpp.html) or downstream
R packages declaring `LinkingTo: rgraphtool, Rcpp, BH`, together with
built-in reference implementations of \\k\\-core decomposition
(`kcore_decomposition()`) and weighted vertex strength
(`vertex_strength()`).

## Usage

``` r
kcore_decomposition(graph, deg = c("total", "out", "in"), checked = FALSE)

vertex_strength(graph, property, mode = c("total", "out", "in"),
  checked = FALSE)

rgraphtool_include_dir()

rgraphtool_extension_example(name = c("kcore", "weighted_degree"))

compile_cpp_extension(file = NULL, code = NULL, env = parent.frame(),
  rebuild = FALSE, verbose = FALSE, ...)
```

## Arguments

- graph:

  A graph created by
  [`graph_create`](https://diogro.github.io/rgraphtool/reference/graph.md).

- deg:

  Degree definition for `kcore_decomposition()`: `"total"` (default,
  matching graph-tool's `kcore_decomposition()`), or `"out"`/`"in"` for
  directed graphs.

- property:

  Character scalar naming an edge property attached to `graph` via
  [`graph_add_edges`](https://diogro.github.io/rgraphtool/reference/graph.md)
  or
  [`graph_set_edge_property`](https://diogro.github.io/rgraphtool/reference/graph.md).

- mode:

  Edge incidence mode for `vertex_strength()`: `"total"`, `"out"`, or
  `"in"` (undirected graphs support only `"total"`).

- checked:

  Logical; if `TRUE`, perform per-access bounds checking on vertex and
  edge property maps in C++; if `FALSE` (default), validate property map
  lengths once via `get_unchecked()` and use zero-overhead direct
  indexing inside hot loops.

- name:

  Bundled C++ extension example name: `"kcore"` or `"weighted_degree"`.

- file:

  Path to a C++ source file (mutually exclusive with `code`).

- code:

  Character vector of inline C++ source code (mutually exclusive with
  `file`).

- env:

  Target R environment where exported C++ functions are assigned.

- rebuild:

  Logical; force recompilation in
  [`sourceCpp`](https://rdrr.io/pkg/Rcpp/man/sourceCpp.html).

- verbose:

  Logical; print compiler invocation details.

- ...:

  Additional arguments forwarded to
  [`sourceCpp`](https://rdrr.io/pkg/Rcpp/man/sourceCpp.html).

## Details

**Public C++ API (`#include <rgraphtool.h>`, namespace `rgraphtool`):**

- *Zero-copy graph access*: `get_graph(SEXP obj)` accepts either an
  `rgraphtool_graph` R object or its external pointer, returning
  `Graph&`. Inspection helpers include `num_vertices(g)`,
  `num_edges(g)`, `is_directed(g)`, `source(e, g)`, `target(e, g)`,
  `edge_index(e)`, `out_degree(v, g)`, `in_degree(v, g)`,
  `total_degree(v, g)`, and `checked_vertex(v, g)`.

- *Range-based graph iteration*: `vertices_range(g)`, `edges_range(g)`,
  `out_edges_range(v, g)`, `in_edges_range(v, g)`,
  `all_edges_range(v, g)`, `out_neighbors_range(v, g)`,
  `in_neighbors_range(v, g)`, and `all_neighbors_range(v, g)`. For
  undirected graphs, `out_edges_range(v, g)` iterates all incident edges
  oriented with `source == v`.

- *Property maps (checked vs. unchecked)*:
  `vertex_property_view(vec, g)` wraps an R `NumericVector` or
  `IntegerVector` zero-copy in a bounds-checked `VertexPropertyMap<T>`;
  `edge_property_view(g, name)` wraps an attached edge property in a
  bounds-checked `EdgePropertyMap<double>` /
  `ConstEdgePropertyMap<double>`; `make_vertex_property<T>(g, init)` and
  `make_edge_property<T>(g, init)` allocate owning property maps.
  Calling `.get_unchecked()` verifies buffer capacity once against
  `num_vertices(g)` / `num_edges(g)` and returns a raw-pointer unchecked
  accessor for tight loops (analogous to graph-tool's
  `PropertyMap.get_unchecked()`).

**Build requirements:** Compiling C++ extensions requires a
C++23-capable toolchain (matching rgraphtool's
`SystemRequirements: C++23`) and the Rcpp and BH header packages.
Downstream R packages can add `LinkingTo: rgraphtool, Rcpp, BH` and
`CXX_STD = CXX23` in `src/Makevars`.

## Value

`kcore_decomposition()` returns an integer vector of length \\N\\
containing each vertex's core number. `vertex_strength()` returns a
numeric vector of length \\N\\. `rgraphtool_include_dir()` and
`rgraphtool_extension_example()` return normalized file paths.
`compile_cpp_extension()` returns the result of
[`sourceCpp`](https://rdrr.io/pkg/Rcpp/man/sourceCpp.html) invisibly.

## Examples

``` r
g <- graph_create(6L, directed = FALSE)
graph_add_edges(g, rbind(c(1, 2), c(2, 3), c(3, 1), c(3, 4), c(4, 5)),
                properties = list(w = c(1, 2, 3, 4, 5)))
kcore_decomposition(g)
#> [1] 2 2 2 1 1 0
vertex_strength(g, "w")
#> [1] 4 3 9 9 5 0
rgraphtool_include_dir()
#> [1] "/private/var/folders/ph/fmqmjzln2h15bm5bstsgpng80000gn/T/Rtmpv0qKtT/temp_libpath170fa7373304a/rgraphtool/include"
rgraphtool_extension_example("kcore")
#> [1] "/private/var/folders/ph/fmqmjzln2h15bm5bstsgpng80000gn/T/Rtmpv0qKtT/temp_libpath170fa7373304a/rgraphtool/examples/extensions/kcore_extension.cpp"
```
