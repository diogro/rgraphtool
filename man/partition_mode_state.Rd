\name{partition_mode_state}
\alias{partition_mode_state}
\alias{mode_marginals}
\alias{mode_consensus}
\alias{mode_uncertainty}
\alias{mode_entropy}
\alias{mode_cluster_state}
\alias{mode_summary}
\alias{align_partition_labels}
\alias{partition_overlap}
\alias{align_nested_partition_labels}
\alias{nested_partition_overlap}
\alias{print.graphtool_partition_mode}
\alias{print.graphtool_mode_cluster}
\title{Partition label alignment, vertex marginals, consensus, and posterior modes}
\usage{
align_partition_labels(x, y)
partition_overlap(x, y, norm = TRUE)
align_nested_partition_labels(x, y)
nested_partition_overlap(x, y, norm = TRUE)
partition_mode_state(partitions, relabel = TRUE, nested = NULL,
  converge = TRUE)
mode_marginals(mode, probability = FALSE, level = 0L)
mode_consensus(mode)
mode_uncertainty(mode, mle = TRUE)
mode_entropy(mode, mle = TRUE)
mode_cluster_state(partitions, nested = NULL, b = NULL, relabel = TRUE,
  niter = 20L, converge = TRUE)
mode_summary(cluster)
}
\arguments{
  \item{x, y}{Flat integer label vectors (1-based) or nested lists of label
  vectors.}
  \item{norm}{Whether to return the normalized overlap in \eqn{[0, 1]} or the raw
  number of matching nodes.}
  \item{partitions}{A list of flat or nested partitions, or the result object
  returned by \code{\link{sample_posterior}}.}
  \item{relabel}{Whether to align partition labels via maximum-weight bipartite
  matching when adding partitions.}
  \item{nested}{Logical indicating whether partitions are hierarchical, or
  \code{NULL} to infer from \code{partitions}.}
  \item{converge}{Whether to iterate label realignment until the random label
  model description length converges.}
  \item{mode}{A partition mode object returned by \code{partition_mode_state()}
  or stored inside \code{mode_cluster_state()$modes}.}
  \item{probability}{If \code{TRUE}, return row-normalized marginal probabilities;
  otherwise return integer counts.}
  \item{level}{Hierarchy level (0-based) for nested modes.}
  \item{mle}{If \code{TRUE}, use maximum-likelihood marginal probabilities; if
  \code{FALSE}, use Dirichlet(1) posterior mean estimates.}
  \item{b}{Optional initial cluster assignment (1-based integer vector of length
  \code{length(partitions)}).}
  \item{niter}{Maximum number of split/sweep/merge optimization iterations. Set
  \code{niter = 0} with \code{b} to score a fixed clustering.}
  \item{cluster}{A mode cluster object returned by \code{mode_cluster_state()}.}
}
\description{
Label alignment, vertex marginal distributions, MAP consensus partitions,
posterior uncertainty, and multimodal clustering of posterior partition ensembles
using Peixoto's random label model (2021). These are the R counterparts of
graph-tool's \code{align_partition_labels()}, \code{partition_overlap()},
\code{PartitionModeState}, and \code{ModeClusterState}.
}
\details{
Pairwise label alignment solves the exact maximum-weight bipartite matching
problem on the contingency matrix using the Hungarian algorithm. For ensembles,
\code{partition_mode_state()} aligns all partitions to a common labelling by
minimizing the integrated description length of the random label model, matching
graph-tool 3.9's \code{PartitionModeState} entropy, \code{posterior_entropy()},
and \code{posterior_cdev()} numerically.

When the posterior distribution is multimodal, \code{mode_cluster_state()}
clusters the sampled partitions into distinct modes by minimizing the mixture
random label model description length, sorts the discovered modes by decreasing
weight, and aligns secondary modes to the dominant mode.
}
\value{
\code{partition_mode_state()} returns a \code{graphtool_partition_mode} list
with \code{M}, \code{B}, \code{entropy}, \code{posterior_entropy},
\code{posterior_cdev}, \code{marginals}, \code{consensus}, and aligned
\code{partitions}. \code{mode_cluster_state()} returns a
\code{graphtool_mode_cluster} list with total \code{entropy}, \code{clusters},
number of modes \code{K}, effective number of modes \code{Be}, \code{weights},
and \code{modes}. \code{mode_summary()} formats the modes as a data frame.
}
\references{
Peixoto, T. P. (2021). Revealing consensus and dissensus between network
partitions. \emph{Physical Review X} 11, 021003.
doi:10.1103/PhysRevX.11.021003
}
\seealso{\code{\link{sample_posterior}}, \code{\link{plot_marginals}},
\code{demo("posterior_sampling")}, \code{demo("hierarchical_sampling")}}
\examples{
bs <- list(c(1, 1, 1, 2, 2, 2), c(2, 2, 2, 1, 1, 1), c(1, 1, 2, 2, 2, 3))
pm <- partition_mode_state(bs)
pm
mode_marginals(pm, probability = TRUE)
mode_consensus(pm)
mode_uncertainty(pm)
}
