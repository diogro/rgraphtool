\name{block_state}
\alias{block_state}
\alias{block_membership}
\alias{block_sizes}
\alias{block_edge_counts}
\alias{block_log_likelihood}
\alias{block_description_length}
\title{Inspect a fixed block-model partition}
\usage{
block_state(graph, blocks = NULL, degree_corrected = TRUE)
block_membership(state)
block_sizes(state)
block_edge_counts(state)
block_log_likelihood(state)
block_description_length(state, degree_prior = NULL)
}
\arguments{
  \item{graph}{A graph created by \code{graph_create()}.}
  \item{blocks}{One positive integer group label per vertex; defaults to a single group.}
  \item{degree_corrected}{Whether to use the degree-corrected likelihood.}
  \item{state}{A state created by \code{block_state()}.}
  \item{degree_prior}{Degree-sequence prior, either \code{"distributed"} or
  \code{"uniform"}. For \code{block_description_length()}, \code{NULL} uses
  the prior associated with the state.}
}
\description{Build and inspect a fixed block-model partition and score it with
the microcanonical stochastic block model description length of graph-tool 3.9.
See \code{\link{minimize_blockmodel_dl}} to fit a partition.}
\details{
Matrix rows and columns correspond to the sorted nonempty labels. Directed
counts are source by target. For undirected graphs the matrix is symmetric,
and diagonal values count each within-group edge twice, including self-loops.
Parallel edges and self-loops are included. The likelihood is the log
probability of the observed multigraph conditional on the block counts,
partition, and, for degree correction, degree sequence. It does not include
priors or description length. \code{block_description_length()} adds the
partition, edge-count, and, for degree-corrected states, degree-sequence
priors, matching graph-tool's \code{BlockState.entropy()} with its default
arguments (\code{dense = FALSE}, \code{multigraph = TRUE}, \code{exact = TRUE}).
The partition prior includes the uniform prior on the number of groups,
\eqn{\log N}.

The distributed degree prior needs the number of integer partitions
\eqn{q(m, n)}. As in graph-tool, exact values are tabulated up to
\eqn{\min(\max(2E, 100), 10000)} and an asymptotic approximation is used
above that. The exact table needs about \eqn{4 n^2} bytes, so this package caps
it at \code{getOption("rgraphtool.q_cache_max", 5000)} (about 100 MB); raise
the option to 10000 for exact agreement with graph-tool whenever a group's
degree sum lies between 5000 and 10000. The table is kept for the session.

States are snapshots: later graph changes do not change their statistics. If
the graph has grown since the state was created, the description length is
computed from the stored snapshot in R (exact partition counts, limited to 10
million table entries) instead of the native code.
}
\value{
\code{block_state()} returns a block state. \code{block_membership()} returns
vertex labels, \code{block_sizes()} returns counts by label,
\code{block_edge_counts()} returns the affinity matrix, and
\code{block_log_likelihood()} returns a numeric scalar in nats.
\code{block_description_length()} returns the description length in nats.
}
\examples{
g <- graph_create(3, directed = FALSE)
graph_add_edges(g, rbind(c(1, 2), c(2, 3)))
state <- block_state(g, c(1, 1, 2))
block_edge_counts(state)
block_log_likelihood(state)
block_description_length(state)
}
\seealso{\code{\link{minimize_blockmodel_dl}}}
