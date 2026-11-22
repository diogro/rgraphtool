\name{block_state}
\alias{block_state}
\alias{block_membership}
\alias{block_sizes}
\alias{block_edge_counts}
\alias{block_log_likelihood}
\alias{block_description_length}
\alias{minimize_blockmodel_dl}
\alias{compare_blockmodels}
\title{Inspect a fixed block-model partition}
\usage{
block_state(graph, blocks = NULL, degree_corrected = TRUE)
block_membership(state)
block_sizes(state)
block_edge_counts(state)
block_log_likelihood(state)
block_description_length(state, degree_prior = NULL)
minimize_blockmodel_dl(graph, degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"), restarts = 3L,
  max_sweeps = 10L)
compare_blockmodels(graph, restarts = 3L, max_sweeps = 10L,
  degree_prior = c("distributed", "uniform"))
}
\arguments{
  \item{graph}{A graph created by \code{graph_create()}.}
  \item{blocks}{One positive integer group label per vertex; defaults to a single group.}
  \item{degree_corrected}{Whether to use the degree-corrected likelihood.}
  \item{state}{A state created by \code{block_state()}.}
  \item{degree_prior}{Degree-sequence prior, either \code{"distributed"} or
  \code{"uniform"}. For \code{block_description_length()}, \code{NULL} uses
  the prior associated with the state.}
  \item{restarts}{Number of search starts, including the single-group start.}
  \item{max_sweeps}{Maximum number of vertex-and-merge sweeps per start.}
}
\description{Build and inspect a fixed block-model partition, or fit a partition
using a greedy minimum-description-length search.}
\details{
Matrix rows and columns correspond to the sorted nonempty labels. Directed
counts are source by target. For undirected graphs the matrix is symmetric,
and diagonal values count each within-group edge twice, including self-loops.
Parallel edges and self-loops are included. The likelihood is the log
probability of the observed multigraph conditional on the block counts,
partition, and, for degree correction, degree sequence. It does not include
priors or description length. \code{block_description_length()} adds the
partition, edge-count, and, for degree-corrected states, degree-sequence
priors. It uses exact integer-partition counts for the distributed prior; the
table is limited to 10 million entries and an error is raised if that limit is
exceeded. The partition prior includes the uniform prior on the number of
groups, \eqn{\log N}.

\code{minimize_blockmodel_dl()} uses greedy single-vertex moves and group
merges from a single-group initialization and random restarts. It uses R's
random-number generator, so \code{set.seed()} controls reproducibility. This is
not graph-tool's multilevel optimizer and can produce worse local optima; it
does not guarantee a global optimum. Its repeated full score computations and
dense matrices make it unsuitable for large graphs.

States are snapshots: later graph changes do not change their statistics.
}
\value{
\code{block_state()} returns a block state. \code{block_membership()} returns
vertex labels, \code{block_sizes()} returns counts by label,
\code{block_edge_counts()} returns the affinity matrix, and
\code{block_log_likelihood()} returns a numeric scalar in nats.
\code{block_description_length()} returns the description length in nats,
and \code{minimize_blockmodel_dl()} returns the best state found.
\code{compare_blockmodels()} returns both fitted states, a score table, and the
minimum-description-length winner.
}
\examples{
g <- graph_create(3, directed = FALSE)
graph_add_edges(g, rbind(c(1, 2), c(2, 3)))
state <- block_state(g, c(1, 1, 2))
block_edge_counts(state)
block_log_likelihood(state)
block_description_length(state)
set.seed(1)
fit <- minimize_blockmodel_dl(g)
block_membership(fit)
}
