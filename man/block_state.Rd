\name{block_state}
\alias{block_state}
\alias{block_membership}
\alias{block_sizes}
\alias{block_edge_counts}
\alias{block_log_likelihood}
\title{Inspect a fixed block-model partition}
\usage{
block_state(graph, blocks = NULL, degree_corrected = TRUE)
block_membership(state)
block_sizes(state)
block_edge_counts(state)
block_log_likelihood(state)
}
\arguments{
  \item{graph}{A graph created by \code{graph_create()}.}
  \item{blocks}{One positive integer group label per vertex; defaults to a single group.}
  \item{degree_corrected}{Whether to use the degree-corrected likelihood.}
  \item{state}{A state created by \code{block_state()}.}
}
\description{Build a snapshot of a fixed partition and inspect its block statistics
and microcanonical SBM log likelihood. This does not infer a partition.}
\details{
Matrix rows and columns correspond to the sorted nonempty labels. Directed
counts are source by target. For undirected graphs the matrix is symmetric,
and diagonal values count each within-group edge twice, including self-loops.
Parallel edges and self-loops are included. The likelihood is the log
probability of the observed multigraph conditional on the block counts,
partition, and, for degree correction, degree sequence. It does not include
priors or description length and cannot be used as a substitute for
\code{graph-tool}'s model-selection objective. States are snapshots: later graph
changes do not change their statistics.
}
\value{
\code{block_state()} returns a block state. \code{block_membership()} returns
vertex labels, \code{block_sizes()} returns counts by label,
\code{block_edge_counts()} returns the affinity matrix, and
\code{block_log_likelihood()} returns a numeric scalar in nats.
}
\examples{
g <- graph_create(3, directed = FALSE)
graph_add_edges(g, rbind(c(1, 2), c(2, 3)))
state <- block_state(g, c(1, 1, 2))
block_edge_counts(state)
block_log_likelihood(state)
}
