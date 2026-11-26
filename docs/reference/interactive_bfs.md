# Interactive BFS shortest-path tree explorer

R-native counterpart of graph-tool's interactive BFS visualization
(`interactive_bst.py`), highlighting the shortest-path BFS tree up to
depth `max_dist` from any selected vertex.

## Usage

``` r
interactive_bfs(graph, pos = NULL, source = 1L, max_dist = 3L,
  directed = NULL, mode = c("html", "locator", "snapshot"), output = NULL,
  width = 760L, height = 560L, highlight_color = "#CE5C00",
  base_color = "#999999", open = FALSE,
  main = "Interactive BFS Tree Explorer")
```

## Arguments

- graph:

  A graph created by
  [`graph_create`](https://diogro.github.io/rgraphtool/reference/graph.md).

- pos:

  Optional \\(N \times 2\\) vertex layout matrix; computed via
  [`sfdp_layout`](https://diogro.github.io/rgraphtool/reference/sfdp_layout.md)
  if `NULL`.

- source:

  1-based integer vertex index (or vector of vertex indices when
  `mode = "snapshot"`) specifying the initial BFS root(s).

- max_dist:

  Maximum BFS depth to highlight (default `3L`, matching graph-tool's
  `interactive_bst.py`).

- directed:

  Logical or `NULL`; defaults to `isTRUE(graph$directed)`.

- mode:

  Interaction mode: `"html"` writes a standalone zero-dependency HTML5
  Canvas web app with live hover/click BFS tree highlighting and depth
  slider; `"locator"` runs an interactive base-R click loop via
  [`locator`](https://rdrr.io/r/graphics/locator.html); `"snapshot"`
  renders static single- or multi-panel BFS tree figure(s).

- output:

  Optional destination HTML path (`mode = "html"`) or PNG file path
  (`mode = "snapshot"`).

- width, height:

  Canvas or PNG image dimensions in pixels.

- highlight_color:

  Colour for BFS tree edges and reached vertices (default `"#CE5C00"`).

- base_color:

  Colour for non-tree edges and unreached vertices.

- open:

  Logical; if `TRUE` in an interactive R session with `mode = "html"`,
  opens the HTML file via
  [`browseURL`](https://rdrr.io/r/utils/browseURL.html).

- main:

  Title for the viewer or plot.

## Details

**Platform and display requirements:**

- `mode = "html"` works headless during generation and runs
  interactively offline in any modern web browser or RStudio Viewer
  (hovering updates the BFS tree live; clicking pins or unpins the
  current root).

- `mode = "locator"` requires an interactive R session with a
  locator-capable graphics device (Quartz, X11, Windows, or RStudioGD);
  it errors cleanly in non-interactive sessions.

- `mode = "snapshot"` works on any R graphics device or writes a PNG
  file via [`png`](https://rdrr.io/r/grDevices/png.html) when `output`
  is supplied.

## Value

An object of class `"rgraphtool_interactive_bfs"` containing `bfs` (an
[`graph_bfs`](https://diogro.github.io/rgraphtool/reference/sfdp_layout.md)
result), `pos`, `mode`, `output`, `max_dist`, and `directed`.

## Examples

``` r
g <- graph_create(10L, directed = FALSE)
graph_add_edges(g, rbind(cbind(1:9, 2:10), c(1L, 6L)))
set.seed(1)
pos <- sfdp_layout(g, max_iter = 30L)
res <- interactive_bfs(g, pos = pos, source = 1L, max_dist = 3L,
                       mode = "html", output = tempfile(fileext = ".html"))
print(res)
#> rgraphtool interactive BFS visualization
#>   mode: html -> /tmp/RtmpTIoacd/filedc83a21b5f4af.html
#> rgraphtool BFS traversal
#>   source: 1 (undirected, max_dist = 3)
#>   reached: 8 / 10 vertices (max distance = 3, tree edges = 7)
```
