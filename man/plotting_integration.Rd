\name{plotting_integration}
\alias{plotting_integration}
\alias{prop_to_size}
\alias{map_colors}
\alias{plot_colorbar}
\alias{save_graph_plot}
\alias{graph_grob}
\alias{graph_to_data_frames}
\alias{ggplot_graph}
\title{Plotting integration: property scaling, colorbars, vector/raster export, grid, and ggplot2}
\usage{
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
}
\arguments{
  \item{x}{Finite numeric vector of vertex or edge property values.}
  \item{mi, ma}{Minimum and maximum target sizes for \code{prop_to_size()}.}
  \item{power}{Positive power-law exponent applied after normalizing \code{x}
  to \code{[0, 1]}.}
  \item{log}{Logical; if \code{TRUE}, log-transform strictly positive values
  before normalizing.}
  \item{palette}{Name of any HCL colour palette supported by
  \code{\link[grDevices]{hcl.colors}}.}
  \item{n}{Number of discrete palette colours.}
  \item{limits}{Numeric length-2 vector \code{c(min, max)} (or, for
  \code{plot_colorbar()}, an object returned by \code{map_colors()}).}
  \item{alpha}{Opacity in \code{[0, 1]}.}
  \item{rev}{Logical; whether to reverse the colour palette.}
  \item{title}{Optional colorbar title.}
  \item{horizontal}{Logical; draw a horizontal colorbar instead of vertical.}
  \item{n_ticks}{Number of axis ticks on the colorbar.}
  \item{digits}{Significant digits for colorbar tick labels.}
  \item{cex_axis}{Tick label character expansion.}
  \item{...}{Additional arguments forwarded to \code{\link[graphics]{axis}}.}
  \item{expr}{A base-R plotting expression, a \code{grid} grob (such as
  \code{graph_grob()}), or a \code{ggplot} object.}
  \item{output}{Destination file path.}
  \item{width, height}{Output dimensions in inches.}
  \item{dpi}{Resolution in dots per inch for raster formats.}
  \item{fmt}{Output format (\code{"auto"} infers from the file extension).
  Vector formats: \code{"pdf"}, \code{"svg"}, \code{"ps"}, \code{"eps"}.
  Raster formats: \code{"png"}, \code{"jpeg"}, \code{"tiff"}.}
  \item{bg}{Background colour.}
  \item{graph}{A graph created by \code{\link{graph_create}}.}
  \item{pos}{Optional \eqn{N \times 2} vertex layout matrix.}
  \item{vertex_color, vertex_border, vertex_size, vertex_halo,
  vertex_halo_color, vertex_halo_size, vertex_visible, vertex_label,
  vertex_label_cex}{Vertex aesthetics matching \code{\link{plot_graph}}.}
  \item{edge_color, edge_alpha, edge_lwd, edge_curve, edge_visible}{Edge
  aesthetics matching \code{\link{plot_graph}}.}
  \item{vorder, eorder}{Optional z-order numeric vectors.}
  \item{arrows}{Logical or \code{NULL}; defaults to
  \code{isTRUE(graph$directed)}.}
  \item{xlim, ylim, pad}{Viewport limits and padding.}
  \item{axes}{Logical; whether to draw coordinate axes.}
  \item{main, subtitle}{Optional plot title and subtitle.}
  \item{vertex_data, edge_data}{Optional data frames or lists of extra vertex
  or edge columns appended by \code{graph_to_data_frames()}.}
  \item{color_label, size_label}{Legend titles for mapped aesthetics in
  \code{ggplot_graph()}.}
}
\description{
R-native plotting integration equivalents of graph-tool's matplotlib workflow:
property-to-size scaling (\code{prop_to_size()}), HCL colour mapping and
standalone colorbars (\code{map_colors()}, \code{plot_colorbar()}), unified
vector/raster export (\code{save_graph_plot()}), embeddable \pkg{grid} grobs
(\code{graph_grob()}), tidy vertex/edge tables
(\code{graph_to_data_frames()}), and optional \pkg{ggplot2} rendering
(\code{ggplot_graph()}).
}
\details{
Vector export supports \code{pdf}, \code{ps}/\code{eps}, and \code{svg} (using
\pkg{svglite} when installed, which avoids XQuartz/Cairo dependencies on macOS,
falling back to \code{\link[grDevices]{svg}}). Raster export supports
\code{png}, \code{jpeg}, and \code{tiff}. \code{ggplot_graph()} requires the
suggested \pkg{ggplot2} package; all other functions rely only on base R
(\pkg{graphics}, \pkg{grDevices}, and \pkg{grid}).
}
\value{
\code{prop_to_size()} returns a numeric vector; \code{map_colors()} returns a
character vector of class \code{"rgraphtool_color_map"}; \code{plot_colorbar()}
invisibly returns the limits; \code{save_graph_plot()} invisibly returns
\code{output}; \code{graph_grob()} returns a \code{grid} \code{"gTree"} of
class \code{"rgraphtool_grob"}; \code{graph_to_data_frames()} returns a list
with \code{vertices} and \code{edges} data frames; \code{ggplot_graph()}
returns a \code{ggplot} object.
}
\examples{
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
}
