\name{plot_bipartite}
\alias{plot_bipartite}
\title{Visualize a fitted bipartite stochastic block model}
\usage{
plot_bipartite(state, type = c("matrix", "meso", "network"), level = 1L,
  weight = NULL, main = NULL, palette = NULL, block_border = "#222222", ...)
}
\arguments{
  \item{state}{An \code{rgraphtool_bipartite_block_state} or
  \code{rgraphtool_nested_bipartite_block_state}.}
  \item{type}{Visualization type: \code{"matrix"} (block-reordered incidence or
  weight heatmap with block boundaries), \code{"meso"} (two-column mesoscale
  block interaction diagram), or \code{"network"} (two-column bipartite
  node-link diagram colored by inferred block).}
  \item{level}{Hierarchy level (1-based) for nested bipartite states.}
  \item{weight}{Optional edge property name to plot in \code{type = "matrix"}.}
  \item{main}{Optional plot title.}
  \item{palette}{Optional color vector for \code{type = "matrix"}.}
  \item{block_border}{Border color for block separators in
  \code{type = "matrix"}.}
  \item{...}{Additional graphical arguments passed to the underlying base-R
  plotting function.}
}
\description{
Produces base-R visualizations of fitted bipartite stochastic block models:
block-reordered incidence matrices, mesoscale block-interaction diagrams, and
two-column bipartite network layouts.
}
\value{
Invisibly returns the plotted summary data (the reordered matrix and orders for
\code{"matrix"}, the density/count matrices for \code{"meso"}, or the 2D vertex
coordinates for \code{"network"}).
}
\examples{
set.seed(1)
M <- matrix(0L, 10, 10)
M[1:5, 1:5] <- rbinom(25, 1, 0.8)
M[6:10, 6:10] <- rbinom(25, 1, 0.8)
fit <- minimize_bipartite_blockmodel_dl(M, dim_labels = c("host", "parasite"))
plot_bipartite(fit, type = "matrix")
plot_bipartite(fit, type = "meso")
}
