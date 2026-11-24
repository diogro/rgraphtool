# Plotting integration: property scaling, colorbars, vector/raster export, grid, and ggplot2

R-native plotting integration equivalents of graph-tool's matplotlib
workflow: property-to-size scaling (`prop_to_size()`), HCL colour
mapping and standalone colorbars (`map_colors()`, `plot_colorbar()`),
unified vector/raster export (`save_graph_plot()`), embeddable grid
grobs (`graph_grob()`), tidy vertex/edge tables
(`graph_to_data_frames()`), and optional ggplot2 rendering
(`ggplot_graph()`).

## Usage

``` r
prop_to_size(x, mi = 0.6, ma = 2.4, power = 0.5, log = FALSE)

map_colors(x, palette = "Viridis", n = 256L, limits = NULL, alpha = 1,
  rev = FALSE)

plot_colorbar(limits, palette = "Viridis", n = 100L, title = NULL, alpha = 1,
  rev = FALSE, horizontal = FALSE, n_ticks = 5L, digits = 2L,
  cex_axis = 0.85, ...)

save_graph_plot(expr, output, width = 7, height = 5, dpi = 150,
  fmt = c("auto", "pdf", "svg", "png", "ps", "eps", "jpeg", "tiff"),
  bg = "white")

graph_grob(graph, pos = NULL, vertex_color = "#4C72B0",
  vertex_border = "#1F2937", vertex_size = 1.0, vertex_halo = FALSE,
  vertex_halo_color = grDevices::adjustcolor("#E63946", 0.4),
  vertex_halo_size = 1.8, vertex_visible = TRUE, vertex_label = NULL,
  vertex_label_cex = 0.7,
  edge_color = grDevices::adjustcolor("#4B5563", 0.45), edge_lwd = 1.0,
  edge_curve = 0, edge_visible = TRUE, vorder = NULL, eorder = NULL,
  arrows = NULL, xlim = NULL, ylim = NULL, pad = 0.08, bg = NULL,
  axes = FALSE, main = NULL)

graph_to_data_frames(graph, pos = NULL, vertex_data = NULL, edge_data = NULL)

ggplot_graph(graph, pos = NULL, vertex_color = "#4C72B0", vertex_size = 2.6,
  vertex_label = NULL, edge_color = "#64748B", edge_alpha = 0.5,
  edge_lwd = 0.5, edge_curve = 0, arrows = NULL, palette = "Viridis",
  color_label = NULL, size_label = NULL, axes = FALSE, main = NULL,
  subtitle = NULL)
```

## Arguments

- x:

  Finite numeric vector of vertex or edge property values.

- mi, ma:

  Minimum and maximum target sizes for `prop_to_size()`.

- power:

  Positive power-law exponent applied after normalizing `x` to `[0, 1]`.

- log:

  Logical; if `TRUE`, log-transform strictly positive values before
  normalizing.

- palette:

  Name of any HCL colour palette supported by
  [`hcl.colors`](https://rdrr.io/r/grDevices/palettes.html).

- n:

  Number of discrete palette colours.

- limits:

  Numeric length-2 vector `c(min, max)` (or, for `plot_colorbar()`, an
  object returned by `map_colors()`).

- alpha:

  Opacity in `[0, 1]`.

- rev:

  Logical; whether to reverse the colour palette.

- title:

  Optional colorbar title.

- horizontal:

  Logical; draw a horizontal colorbar instead of vertical.

- n_ticks:

  Number of axis ticks on the colorbar.

- digits:

  Significant digits for colorbar tick labels.

- cex_axis:

  Tick label character expansion.

- ...:

  Additional arguments forwarded to
  [`axis`](https://rdrr.io/r/graphics/axis.html).

- expr:

  A base-R plotting expression, a `grid` grob (such as `graph_grob()`),
  or a `ggplot` object.

- output:

  Destination file path.

- width, height:

  Output dimensions in inches.

- dpi:

  Resolution in dots per inch for raster formats.

- fmt:

  Output format (`"auto"` infers from the file extension). Vector
  formats: `"pdf"`, `"svg"`, `"ps"`, `"eps"`. Raster formats: `"png"`,
  `"jpeg"`, `"tiff"`.

- bg:

  Background colour.

- graph:

  A graph created by
  [`graph_create`](https://diogro.github.io/rgraphtool/reference/graph.md).

- pos:

  Optional \\N \times 2\\ vertex layout matrix.

- vertex_color, vertex_border, vertex_size, vertex_halo,
  vertex_halo_color, vertex_halo_size, vertex_visible, vertex_label,
  vertex_label_cex:

  Vertex aesthetics matching
  [`plot_graph`](https://diogro.github.io/rgraphtool/reference/plot_graph.md).

- edge_color, edge_alpha, edge_lwd, edge_curve, edge_visible:

  Edge aesthetics matching
  [`plot_graph`](https://diogro.github.io/rgraphtool/reference/plot_graph.md).

- vorder, eorder:

  Optional z-order numeric vectors.

- arrows:

  Logical or `NULL`; defaults to `isTRUE(graph$directed)`.

- xlim, ylim, pad:

  Viewport limits and padding.

- axes:

  Logical; whether to draw coordinate axes.

- main, subtitle:

  Optional plot title and subtitle.

- vertex_data, edge_data:

  Optional data frames or lists of extra vertex or edge columns appended
  by `graph_to_data_frames()`.

- color_label, size_label:

  Legend titles for mapped aesthetics in `ggplot_graph()`.

## Details

Vector export supports `pdf`, `ps`/`eps`, and `svg` (using svglite when
installed, which avoids XQuartz/Cairo dependencies on macOS, falling
back to [`svg`](https://rdrr.io/r/grDevices/cairo.html)). Raster export
supports `png`, `jpeg`, and `tiff`. `ggplot_graph()` requires the
suggested ggplot2 package; all other functions rely only on base R
(graphics, grDevices, and grid).

## Value

`prop_to_size()` returns a numeric vector; `map_colors()` returns a
character vector of class `"rgraphtool_color_map"`; `plot_colorbar()`
invisibly returns the limits; `save_graph_plot()` invisibly returns
`output`; `graph_grob()` returns a `grid` `"gTree"` of class
`"rgraphtool_grob"`; `graph_to_data_frames()` returns a list with
`vertices` and `edges` data frames; `ggplot_graph()` returns a `ggplot`
object.

## Examples

``` r
g <- graph_create(8L, directed = FALSE)
graph_add_edges(g, rbind(cbind(1:7, 2:8), c(8L, 1L), c(1L, 5L)))
set.seed(1)
pos <- sfdp_layout(g, max_iter = 30L)
deg <- graph_degree(g)
vcol <- map_colors(deg, palette = "Plasma")
vsiz <- prop_to_size(deg, mi = 0.8, ma = 2.0)
gb <- graph_grob(g, pos = pos, vertex_color = vcol, vertex_size = vsiz,
                 axes = TRUE, main = "Grid Graph")
save_graph_plot(gb, output = tempfile(fileext = ".pdf"), width = 5, height = 4)
```
