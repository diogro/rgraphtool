\name{pagerank}
\alias{pagerank}
\alias{betweenness}
\alias{closeness}
\alias{eigenvector_centrality}
\alias{katz_centrality}
\alias{hits_centrality}
\title{Vertex and edge centrality measures}
\usage{
pagerank(graph, damping = 0.85, pers = NULL, weights = NULL,
         tol = 1e-6, max_iter = 1000L)
betweenness(graph, weights = NULL, normalized = TRUE)
closeness(graph, weights = NULL, normalized = TRUE, harmonic = FALSE)
eigenvector_centrality(graph, weights = NULL, tol = 1e-6, max_iter = 1000L)
katz_centrality(graph, alpha = 0.01, beta = NULL, weights = NULL,
                normalized = TRUE, tol = 1e-6, max_iter = 1000L)
hits_centrality(graph, weights = NULL, tol = 1e-6, max_iter = 1000L)
}
\arguments{
  \item{graph}{An \code{rgraphtool_graph} object.}
  \item{damping}{PageRank damping factor in \eqn{[0, 1]} (default \code{0.85}).}
  \item{pers}{Optional non-negative personalization vector of length
  \code{graph_num_vertices(graph)}.}
  \item{weights}{Optional edge property name or numeric vector of length
  \code{graph_num_edges(graph)}. For \code{betweenness()} and
  \code{closeness()}, weights are interpreted as path lengths (distances); for
  spectral/walk centralities, weights are interpreted as connection strengths.}
  \item{tol}{Convergence tolerance for iterative power methods.}
  \item{max_iter}{Maximum number of power iterations.}
  \item{normalized}{Logical; whether to normalize centrality values.}
  \item{harmonic}{Logical; if \code{TRUE}, computes harmonic closeness
  \eqn{\sum_{v \neq u} 1 / d(u, v)} (divided by \eqn{N - 1} when
  \code{normalized = TRUE}).}
  \item{alpha}{Katz attenuation factor (should be smaller than
  \eqn{1 / \lambda_{\max}}).}
  \item{beta}{Optional exogenous vertex weight vector of length \eqn{N} for
  \code{katz_centrality()}.}
}
\description{
Multi-threaded native C++23 implementations of classic network centrality
measures matching \code{graph-tool}'s \code{centrality} module: PageRank,
Brandes vertex and edge betweenness, closeness (standard and harmonic), leading
eigenvector centrality, Katz centrality, and Kleinberg HITS hubs and
authorities.
}
\value{
\code{pagerank()}, \code{closeness()}, and \code{katz_centrality()} return
numeric vectors of length \eqn{N}. \code{betweenness()} returns a list with
elements \code{vertex} (length \eqn{N}) and \code{edge} (length \eqn{E}).
\code{eigenvector_centrality()} returns \code{list(eigenvalue, vector)}, and
\code{hits_centrality()} returns \code{list(eigenvalue, authority, hub)}.
}
\examples{
g <- graph_create(6L, directed = FALSE)
graph_add_edges(g, rbind(c(1, 2), c(2, 3), c(3, 1), c(3, 4), c(4, 5), c(5, 6), c(6, 4)))
pagerank(g)
betweenness(g)
closeness(g)
}
