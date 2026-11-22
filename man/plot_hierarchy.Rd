\name{plot_hierarchy}
\alias{plot_hierarchy}
\title{Draw a nested block model as a radial hierarchy}
\usage{
plot_hierarchy(state, beta = 0.8, vertex_col = NULL,
  edge_col = grDevices::adjustcolor("grey30", 0.25), vertex_cex = 0.6,
  show_tree = TRUE, main = NULL, ...)
}
\arguments{
  \item{state}{A nested block state.}
  \item{beta}{Edge bundling strength between 0 (straight lines) and 1 (edges
  follow the hierarchy).}
  \item{vertex_col}{Vertex colours; by default, one colour per level-0 group.}
  \item{edge_col}{Edge colour, recycled over edges.}
  \item{vertex_cex}{Vertex symbol size.}
  \item{show_tree}{Whether to draw the hierarchy tree.}
  \item{main}{Optional title.}
  \item{...}{Passed to \code{points()} for the vertices.}
}
\description{
A base-R counterpart of graph-tool's \code{NestedBlockState.draw()}. Vertices lie
on the unit circle, ordered by the hierarchy. Groups sit at radii that shrink
towards the root at the centre. Edges are drawn as splines through their
endpoints' common ancestors (hierarchical edge bundling, Holten 2006).
}
\details{
The layout differs from graph-tool's \code{draw_hierarchy()}. Positions are
deterministic and evenly spaced rather than computed with its radial tree
layout, and edge direction is not shown. Self-loops are omitted.
}
\value{Invisibly, the layout: vertex coordinates, group coordinates by level,
and each vertex's ancestors.}
\examples{
g <- graph_create(12, directed = FALSE)
graph_add_edges(g, do.call(rbind, lapply(list(1:4, 5:8, 9:12),
                                         function(v) t(combn(v, 2)))))
state <- nested_block_state(g, list(rep(1:3, each = 4), c(1, 1, 2), c(1, 1)))
plot_hierarchy(state)
}
