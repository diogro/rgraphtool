\name{plot_marginals}
\alias{plot_marginals}
\title{Draw a network with pie-chart vertex marginal group probabilities}
\usage{
plot_marginals(graph, mode, positions = NULL, radius = 0.04,
  palette = NULL, edge_col = grDevices::adjustcolor("grey40", 0.35),
  border = "grey20", labels = NULL, label_cex = 0.65, main = NULL, ...)
}
\arguments{
  \item{graph}{A graph created by \code{\link{graph_create}}.}
  \item{mode}{A \code{graphtool_partition_mode} from
  \code{\link{partition_mode_state}} or \code{\link{mode_cluster_state}}, or an
  \eqn{N \times B} numeric matrix of marginal counts or probabilities.}
  \item{positions}{Optional \eqn{N \times 2} matrix of vertex coordinates. When
  \code{NULL}, computes a deterministic spring layout initialized from the
  consensus partition.}
  \item{radius}{Pie-chart radius in normalized plot coordinates.}
  \item{palette}{Optional vector of group colours.}
  \item{edge_col}{Edge colour.}
  \item{border}{Pie slice border colour.}
  \item{labels}{Optional character vector of vertex labels.}
  \item{label_cex}{Character expansion for \code{labels}.}
  \item{main}{Optional plot title.}
  \item{...}{Additional arguments passed to \code{graphics::plot.window()}.}
}
\description{
Base-R visualization of posterior vertex marginal group memberships, drawing a
pie chart at each vertex proportional to its marginal probabilities across
aligned groups (the R equivalent of graph-tool's \code{vertex_shape="pie"}
marginal plots).
}
\value{Invisibly, a list with normalized vertex \code{positions}, the marginal
\code{probabilities} matrix, and the group colour \code{palette}.}
\seealso{\code{\link{partition_mode_state}}, \code{\link{plot_hierarchy}},
\code{demo("posterior_sampling")}}
\examples{
g <- graph_create(6, directed = FALSE)
graph_add_edges(g, rbind(c(1, 2), c(2, 3), c(1, 3), c(4, 5), c(5, 6), c(4, 6), c(3, 4)))
pm <- partition_mode_state(list(c(1, 1, 1, 2, 2, 2), c(1, 1, 2, 2, 2, 2)))
plot_marginals(g, pm, main = "Vertex marginals")
}
