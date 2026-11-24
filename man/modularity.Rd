\name{modularity}
\alias{modularity}
\alias{maximize_modularity}
\alias{print.rgraphtool_modularity_state}
\title{Newman modularity and modularity maximization (use with extreme caution)}
\usage{
modularity(graph, blocks, gamma = 1.0, weights = NULL)
maximize_modularity(graph, gamma = 1.0, weights = NULL, nruns = 1L,
                    max_levels = 10L, warn = TRUE)
\method{print}{rgraphtool_modularity_state}(x, ...)
}
\arguments{
  \item{graph}{An \code{rgraphtool_graph} object.}
  \item{blocks}{Integer partition vector of length
  \code{graph_num_vertices(graph)} (or a fitted block/modularity state).}
  \item{gamma}{Resolution parameter \eqn{\gamma} (default \code{1.0}).}
  \item{weights}{Optional edge property name or numeric edge weight vector.}
  \item{nruns}{Number of randomized multilevel Louvain restarts; keeps the
  highest modularity score.}
  \item{max_levels}{Maximum number of Louvain coarsening passes.}
  \item{warn}{Logical; if \code{TRUE} (the default), emits an explicit warning
  explaining why modularity maximization is statistically unsound.}
  \item{x}{An \code{rgraphtool_modularity_state} object.}
  \item{...}{Unused; included for S3 \code{print} compatibility.}
}
\description{
Computes generalized Newman--Girvan / Reichardt--Bornholdt modularity
\eqn{Q(\gamma)} and maximizes it via fast native C++23 multilevel Louvain
coarsening and fine-grained refinement.
}
\section{DANGER: Modularity Maximization Considered Harmful}{
Using modularity maximization for community detection is almost always
\strong{a terrible idea}. Following \code{graph-tool}'s upstream documentation,
\code{maximize_modularity()} is included \emph{strictly for comparison and
pedagogical purposes} and emits an explicit warning by default:
\itemize{
  \item \strong{Massive overfitting of random noise}: Modularity has no
  statistical regularization or prior penalty on the number of groups. On
  completely structureless Erdos--Renyi or configuration-model random graphs,
  \code{maximize_modularity()} routinely fabricates spurious partitions with
  high modularity (\eqn{Q > 0.3}--\eqn{0.8}), whereas
  \code{\link{minimize_blockmodel_dl}()} rightly infers \eqn{B = 1}.
  \item \strong{Resolution limit (underfitting)}: Simultaneously, modularity
  cannot resolve genuine communities smaller than \eqn{\sim \sqrt{2E}} edges
  (Fortunato & Barthelemy 2007), merging distinct small modules in large
  networks.
  \item \strong{Inferential alternative}: For principled Bayesian community
  detection—including specifically assortative communities—use
  \code{\link{minimize_blockmodel_dl}()},
  \code{\link{minimize_nested_blockmodel_dl}()}, or
  \code{\link{minimize_pp_blockmodel_dl}()} instead.
}
See Tiago P. Peixoto (2023), \emph{Descriptive vs. inferential community
detection in networks: pitfalls, myths and half-truths} (Cambridge University
Press) and \url{https://skewed.de/tiago/blog/modularity-harmful}.
}
\value{
\code{modularity()} returns the scalar modularity \eqn{Q(\gamma)}.
\code{maximize_modularity()} returns an \code{rgraphtool_modularity_state} list
containing \code{blocks}, \code{modularity}, \code{B}, \code{gamma}, and
\code{graph}.
}
\examples{
# Illustration of modularity overfitting a completely random Erdos-Renyi graph:
set.seed(42)
g_rand <- random_graph_erdos_renyi(80L, p = 0.06)

# Modularity maximization fabricates spurious groups (and warns!):
mod_fit <- suppressWarnings(maximize_modularity(g_rand))
mod_fit$B
mod_fit$modularity

# Principled MDL inference rightly finds B = 1 (no community structure):
sbm_fit <- minimize_blockmodel_dl(g_rand)
length(block_sizes(sbm_fit))
}
