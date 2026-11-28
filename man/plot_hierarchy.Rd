\name{plot_hierarchy}
\alias{plot_hierarchy}
\alias{plot_ranked_hierarchy}
\title{Draw a nested or ranked block model hierarchy}
\usage{
plot_hierarchy(state, beta = 0.8, vertex_col = NULL,
  edge_col = grDevices::adjustcolor("grey30", 0.25), vertex_cex = 0.6,
  show_tree = TRUE, main = NULL, tree_col = "#729fcf", edge_lwd = 1,
  arc_level = NULL, arc_col = NULL, arc_labels = TRUE, arc_radius = 1.06,
  arc_lwd = 5, arc_label_cex = 0.8, arc_label_radius = arc_radius + 0.05, ...)
plot_ranked_hierarchy(state,
  upstream_col = grDevices::adjustcolor("#4878CF", 0.35),
  downstream_col = grDevices::adjustcolor("#D65F5F", 0.55),
  lateral_col = grDevices::adjustcolor("#333355", 0.25),
  vertex_col = NULL, vertex_cex = 0.85, curve = 0.15,
  show_ranks = TRUE, show_legend = TRUE, main = NULL, ...)
}
\arguments{
  \item{state}{For \code{plot_hierarchy()}, a nested block state (including
  nested weighted, layered, and ranked states). For
  \code{plot_ranked_hierarchy()}, a flat or nested ranked block state from
  \code{\link{ranked_block_state}}, \code{\link{minimize_ranked_blockmodel_dl}},
  \code{\link{nested_ranked_block_state}}, or
  \code{\link{minimize_nested_ranked_blockmodel_dl}}.}
  \item{beta}{Edge bundling strength between 0 (straight lines) and 1 (edges
  follow the hierarchy).}
  \item{vertex_col}{Vertex colours; by default, one colour per level-0 group.}
  \item{edge_col}{Edge colour, recycled over edges.}
  \item{vertex_cex}{Vertex symbol size.}
  \item{show_tree}{Whether to draw the hierarchy tree.}
  \item{edge_lwd}{Edge line width, recycled over edges like \code{edge_col}.
    Thin edges are drawn first so heavy ones stay on top.}
  \item{arc_level}{Hierarchy level whose groups are marked by arcs outside the
    circle: \code{0} for the groups of the vertices, \code{1} for the next level
    up, and so on. \code{NULL} (default) draws no arcs.}
  \item{arc_col}{Arc colours: one per group (named by group id, or in group
    order), recycled. By default level-0 arcs repeat the default vertex
    colours and coarser levels alternate two greys.}
  \item{arc_labels}{\code{TRUE} labels each arc with its group id,
    \code{FALSE} draws no labels, and a vector with one value per vertex
    labels each arc with the value at its vertices.}
  \item{arc_radius, arc_lwd}{Radius and line width of the arcs.}
  \item{arc_label_cex, arc_label_radius}{Size and radius of the arc labels.}
  \item{tree_col}{Colour of the hierarchy tree edges and group nodes, drawn at
  alpha 0.6 (default \code{"#729fcf"}, the light blue of graph-tool's
  \code{draw_hierarchy()}).}
  \item{upstream_col, downstream_col, lateral_col}{Colours for upstream
  (\eqn{u_i < u_j}), downstream (\eqn{u_i > u_j}), and lateral
  (\eqn{u_i = u_j}) directed edges in \code{plot_ranked_hierarchy()}.}
  \item{curve}{Curvature applied to directed edge splines in
  \code{plot_ranked_hierarchy()} (\code{0} draws straight segments).}
  \item{show_ranks}{Whether to draw horizontal guide lines at each inferred
  rank tier.}
  \item{show_legend}{Whether to display a legend with upstream, lateral, and
  downstream edge counts.}
  \item{main}{Optional title.}
  \item{...}{Passed to \code{points()} for the vertices.}
}
\description{
\code{plot_hierarchy()} is a base-R counterpart of graph-tool's
\code{NestedBlockState.draw()}. Vertices lie on the unit circle, ordered by the
hierarchy. Groups sit at radii that shrink towards the root at the centre. Edges
are drawn as splines through their endpoints' common ancestors (hierarchical
edge bundling, Holten 2006).

\code{plot_ranked_hierarchy()} draws a vertical rank-ordered layout for flat and
nested ranked stochastic block models, placing vertices at vertical positions
determined by their inferred group rank \eqn{u_r} and colouring curved directed
arcs by flow direction (upstream, lateral, and downstream).
}
\details{
In \code{plot_hierarchy()}, vertex positions are deterministic and evenly spaced
rather than computed with graph-tool's radial tree layout, and edge direction is
not shown. In both layouts, self-loops are omitted.
}
\value{Invisibly, the layout coordinates and group metadata. With
\code{arc_level}, \code{$arcs} lists each group's start, end and middle
angle, size, colour and label.}
\examples{
g <- graph_create(12, directed = FALSE)
graph_add_edges(g, do.call(rbind, lapply(list(1:4, 5:8, 9:12),
                                         function(v) t(combn(v, 2)))))
state <- nested_block_state(g, list(rep(1:3, each = 4), c(1, 1, 2), c(1, 1)))
plot_hierarchy(state)
plot_hierarchy(state, arc_level = 1)
}
