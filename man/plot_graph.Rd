\name{plot_graph}
\alias{plot_graph}
\title{Draw a graph with base R graphics}
\usage{
plot_graph(graph, pos = NULL, vertex_color = "#4C72B0",
  vertex_border = "#1F2937", vertex_size = 1.2, vertex_shape = 21L,
  vertex_halo = FALSE,
  vertex_halo_color = grDevices::adjustcolor("#E63946", 0.45),
  vertex_halo_size = 2.2, vertex_visible = TRUE, vertex_label = NULL,
  vertex_label_cex = 0.7,
  edge_color = grDevices::adjustcolor("#4B5563", 0.45), edge_lwd = 1.0,
  edge_curve = 0, edge_visible = TRUE, vorder = NULL, eorder = NULL,
  arrows = NULL, arrow_length = 0.08, xlim = NULL, ylim = NULL, pad = 0.08,
  bg = "white", axes = FALSE, frame.plot = axes, grid_lines = FALSE,
  xlab = "", ylab = "", main = NULL, add = FALSE, output = NULL, width = 7,
  height = 5, dpi = 150, ...)
}
\arguments{
  \item{graph}{A graph created by \code{\link{graph_create}}.}
  \item{pos}{Optional \eqn{N \times 2} coordinate matrix; computed via
  \code{\link{sfdp_layout}} if \code{NULL}.}
  \item{vertex_color}{Vertex fill colour(s), recycled to length \eqn{N}.}
  \item{vertex_border}{Vertex border colour(s) for filled symbols
  (\code{pch \%in\% 21:25}).}
  \item{vertex_size}{Vertex symbol size(s) (\code{cex}).}
  \item{vertex_shape}{Vertex plotting symbol(s) (\code{pch}).}
  \item{vertex_halo}{Logical vector indicating which vertices receive a
  highlight halo behind the vertex marker.}
  \item{vertex_halo_color}{Halo colour(s).}
  \item{vertex_halo_size}{Halo symbol size(s).}
  \item{vertex_visible}{Logical mask hiding vertices (and their incident edges)
  when \code{FALSE}.}
  \item{vertex_label}{Optional character vector of vertex labels.}
  \item{vertex_label_cex}{Character expansion factor for \code{vertex_label}.}
  \item{edge_color}{Edge colour(s), recycled to \code{graph_num_edges(graph)}.}
  \item{edge_lwd}{Edge line width(s).}
  \item{edge_curve}{Numeric edge curvature (\code{0} draws straight segments;
  non-zero draws quadratic Bezier curves).}
  \item{edge_visible}{Logical mask hiding edges when \code{FALSE}.}
  \item{vorder, eorder}{Optional numeric ordering vectors controlling vertex
  and edge z-order (higher values are drawn on top).}
  \item{arrows}{Logical or \code{NULL}; whether to draw directed arrowheads.
  Defaults to \code{isTRUE(graph$directed)}.}
  \item{arrow_length}{Arrowhead length in inches passed to
  \code{\link[graphics]{arrows}}.}
  \item{xlim, ylim}{Optional axis limits.}
  \item{pad}{Fractional viewport padding added around \code{pos} when
  \code{xlim} or \code{ylim} is \code{NULL}.}
  \item{bg}{Background colour.}
  \item{axes}{Logical; whether to draw coordinate axes.}
  \item{frame.plot}{Logical; whether to draw a box around the plot region.}
  \item{grid_lines}{Logical; whether to draw faint coordinate grid lines.}
  \item{xlab, ylab}{Optional axis labels when \code{axes = TRUE}.}
  \item{main}{Optional plot title.}
  \item{add}{Logical; if \code{TRUE}, draw on top of the current plot.}
  \item{output}{Optional output file path; when non-\code{NULL}, writes
  directly via \code{\link{save_graph_plot}} (PDF, SVG, EPS, PNG, JPEG, TIFF).}
  \item{width, height}{Output figure dimensions in inches when \code{output}
  is supplied.}
  \item{dpi}{Raster resolution in dots per inch when \code{output} is a raster
  format.}
  \item{...}{Additional arguments forwarded to \code{\link[graphics]{points}}.}
}
\description{
Draws a graph using base R graphics without external dependencies, supporting
vertex and edge filtering, highlight halos, curved edges, z-ordering, axes,
directed arrows, and direct vector or raster file export.
}
\value{Invisibly, the \eqn{N \times 2} vertex coordinate matrix used for
plotting.}
\examples{
g <- graph_create(6L, directed = FALSE)
graph_add_edges(g, cbind(c(1, 2, 3, 4, 5, 6), c(2, 3, 4, 5, 6, 1)))
set.seed(1)
plot_graph(g, vertex_halo = c(TRUE, rep(FALSE, 5L)), edge_curve = 0.1,
           axes = TRUE, main = "Hexagon")
}
