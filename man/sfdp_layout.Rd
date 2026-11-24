\name{sfdp_layout}
\alias{sfdp_layout}
\alias{graph_bfs}
\title{Force-directed SFDP layout and breadth-first search traversal}
\usage{
sfdp_layout(graph, pos = NULL, C = 0.2, K = 1.0, p = 2.0, theta = 0.6,
  gamma = 1.0, mu = 0.0, mu_p = 1.0, init_step = NULL, cooling_step = 0.95,
  max_level = 10L, tol = 1e-4, max_iter = 200L, simple = FALSE,
  groups = NULL, pin = NULL)

graph_bfs(graph, source, max_dist = Inf, directed = NULL)
}
\arguments{
  \item{graph}{A graph created by \code{\link{graph_create}}.}
  \item{pos}{Optional \eqn{N \times 2} numeric matrix of initial vertex
  coordinates. Useful for warm-starting or incremental layout updates
  (\code{max_iter = 1L}).}
  \item{C}{Relative repulsive force strength (positive scalar).}
  \item{K}{Optimal spring length (positive scalar).}
  \item{p}{Repulsive force exponent (positive scalar; default \code{2.0}).}
  \item{theta}{Barnes-Hut quadtree opening angle parameter (used when
  \eqn{N > 128} and \code{max_level > 0}).}
  \item{gamma}{Inter-group force coefficient (\code{> 0} attracts distinct
  groups toward each other, \code{< 0} repels them).}
  \item{mu}{Intra-group attractive force coefficient (\code{>= 0}).}
  \item{mu_p}{Intra-group attractive scaling exponent.}
  \item{init_step}{Optional positive initial displacement step size. Defaults
  to \code{max(K, 0.1 * K * sqrt(N))}.}
  \item{cooling_step}{Multiplicative step-size cooling factor in \code{(0, 1]}.}
  \item{max_level}{Maximum Barnes-Hut quadtree depth (integer \code{>= 0}).}
  \item{tol}{Stopping tolerance on mean vertex displacement relative to
  \code{K}.}
  \item{max_iter}{Maximum number of layout iterations (\code{0L} means no cap).}
  \item{simple}{Logical; if \code{TRUE}, use monotone cooling
  (\code{step <- step * cooling_step}); if \code{FALSE}, use Hu (2005) adaptive
  step sizing.}
  \item{groups}{Optional vector of length \eqn{N} assigning vertices to groups
  for inter- and intra-group layout forces.}
  \item{pin}{Optional logical vector of length \eqn{N}; vertices with
  \code{TRUE} remain fixed at their initial \code{pos} coordinates.}
  \item{source}{1-based integer index of the starting vertex for BFS.}
  \item{max_dist}{Maximum BFS depth (\code{Inf} for unrestricted traversal, or
  a non-negative integer).}
  \item{directed}{Logical or \code{NULL}. If \code{NULL}, defaults to
  \code{isTRUE(graph$directed)}.}
}
\description{
\code{sfdp_layout()} computes a two-dimensional spring-electrical
force-directed layout in native C++ following Hu (2005) and graph-tool's
\code{sfdp_layout()}, using exact pairwise repulsion for small graphs
(\eqn{N \le 128}) and a Barnes-Hut quadtree for larger graphs.

\code{graph_bfs()} computes a breadth-first search traversal and shortest-path
spanning tree from \code{source} up to depth \code{max_dist}.
}
\details{
Both functions are deterministic given R's random number seed
(\code{\link[base]{set.seed}} for \code{sfdp_layout()}, and completely
deterministic for \code{graph_bfs()}). Self-loops are ignored during layout and
BFS traversal; parallel edges are collapsed to a single tree edge in
\code{graph_bfs()}.
}
\value{
\code{sfdp_layout()} returns an \eqn{N \times 2} numeric matrix with column
names \code{"x"} and \code{"y"}.

\code{graph_bfs()} returns an object of class \code{"rgraphtool_bfs"} with
elements \code{source}, \code{max_dist}, \code{directed}, \code{distance}
(numeric vector of length \eqn{N}, \code{Inf} for unreached vertices),
\code{parent} (integer vector of 1-based predecessors, \code{NA_integer_} for
the source and unreached vertices), \code{order} (1-based discovery order), and
\code{tree_edges} (logical vector aligned with \code{\link{graph_edges}}).
}
\examples{
g <- graph_create(10L, directed = FALSE)
graph_add_edges(g, cbind(1:9, 2:10))
set.seed(1)
pos <- sfdp_layout(g, max_iter = 40L)
bfs <- graph_bfs(g, source = 1L, max_dist = 3L)
print(bfs)
}
