\name{weighted_block_state}
\alias{weighted_block_state}
\alias{weighted_description_length}
\alias{minimize_weighted_blockmodel_dl}
\alias{nested_weighted_block_state}
\alias{nested_weighted_level_entropy}
\alias{nested_weighted_description_length}
\alias{minimize_nested_weighted_blockmodel_dl}
\alias{print.rgraphtool_weighted_block_state}
\alias{print.rgraphtool_nested_weighted_block_state}
\title{Flat and nested weighted stochastic block models with edge covariates}
\usage{
weighted_block_state(graph, weights, types, params = "microcanonical",
  blocks = NULL, degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"))
weighted_description_length(state, rec = TRUE, rec_dl = TRUE,
  degree_prior = NULL)
minimize_weighted_blockmodel_dl(graph, weights, types,
  params = "microcanonical", degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"), B_min = 1L, B_max = Inf,
  refine = TRUE, epsilon = 1e-3, nthreads = NULL, spectral_starts = TRUE)
nested_weighted_block_state(graph, weights, types, params = "microcanonical",
  levels = NULL, degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"))
nested_weighted_level_entropy(state, degree_prior = NULL)
nested_weighted_description_length(state, degree_prior = NULL)
minimize_nested_weighted_blockmodel_dl(graph, weights, types,
  params = "microcanonical", degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"), B_min = 1L, B_max = Inf,
  refine = TRUE, epsilon = 1e-3, nthreads = NULL, spectral_starts = TRUE)
}
\arguments{
  \item{graph}{A graph created by \code{\link{graph_create}}.}
  \item{weights}{Edge covariate values: either a character vector naming one or
  more edge properties stored on \code{graph}, a numeric vector of length
  \code{graph_num_edges(graph)}, an \eqn{E \times K} numeric matrix, or a list
  mixing property names and numeric vectors.}
  \item{types}{Character vector of length \eqn{K} naming the covariate model for
  each column of \code{weights}: \code{"real-exponential"},
  \code{"real-normal"}, \code{"discrete-geometric"}, \code{"discrete-poisson"},
  \code{"discrete-binomial"}, or \code{"discrete-negbinomial"} (overdispersed
  counts; see Details).}
  \item{params}{Covariate hyperprior specification: \code{"microcanonical"}
  (default, Peixoto 2018 maximum-entropy microcanonical priors),
  \code{"bayesian"} (empirical-Bayes conjugate hyperpriors matching graph-tool),
  a named list of hyperparameters (e.g. \code{list(alpha = 1, beta = 2)} or
  \code{list(N = 10, alpha = 1, beta = 1)} or
  \code{list(m0 = 0, k0 = 1, v0 = 1, nu0 = 3)}; \code{list(k = 2)} or
  \code{list(k = 2, alpha = 1, beta = 1)} for \code{"discrete-negbinomial"},
  which has no default size \code{k}), or a list of length \eqn{K}
  specifying one entry per covariate.}
  \item{blocks}{Optional 1-based integer group membership vector of length
  \code{graph_num_vertices(graph)}.}
  \item{levels}{Optional non-empty list of 1-based integer label vectors for
  each hierarchy level.}
  \item{degree_corrected}{Whether to use the degree-corrected (\code{TRUE}) or
  non-degree-corrected (\code{FALSE}) topology model.}
  \item{degree_prior}{Degree sequence hyperprior (\code{"distributed"} or
  \code{"uniform"}).}
  \item{state}{A \code{rgraphtool_weighted_block_state} or
  \code{rgraphtool_nested_weighted_block_state} object.}
  \item{rec}{Whether to include the edge-covariate log-likelihood term.}
  \item{rec_dl}{Whether to include the global covariate hyperprior description
  length term.}
  \item{B_min, B_max}{Minimum and maximum number of occupied level-0 groups.}
  \item{refine}{Whether to run greedy local refinement sweeps.}
  \item{epsilon}{Convergence tolerance in nats.}
  \item{spectral_starts}{Whether to add weight-informed starting partitions
  (default \code{TRUE}). The agglomerative search proposes merges from
  neighbouring groups, which carry no information when only the edge weights
  are structured. Extra starts are built from the leading eigenvectors of the
  weight-centred adjacency (log1p for non-negative types), by varimax rotation
  (Rohe and Zeng 2023) and k-means for 2 to 8 groups. Each is refined, and the
  lowest description length is kept, so results are never worse; R's random
  state is left unchanged. On graphs whose blocks show only in negative-binomial
  weights, NMI with the planted partition rose from 0.16 to 0.71. Runtime
  roughly doubles on large graphs. Uses the first weight covariate.}
  \item{nthreads}{Optional thread count override (defaults to
  \code{rgraphtool_num_threads()}).}
}
\description{
Construct, evaluate, and infer flat and hierarchical (nested) weighted
stochastic block models (Peixoto 2018) supporting all five graph-tool edge
covariate distributions (\code{real-exponential}, \code{real-normal},
\code{discrete-geometric}, \code{discrete-poisson}, and
\code{discrete-binomial}), multiple simultaneous covariates, and both
microcanonical and Bayesian conjugate hyperpriors.
}
\details{
\code{"discrete-negbinomial"} is not in graph-tool. It models counts as Poisson
with an observation-level random effect: each edge's rate is multiplied by an
independent Gamma(\eqn{k}, \eqn{k}) variable, so counts are negative binomial
with variance \eqn{\mu + \mu^2 / k}. Microcanonically, the counts of a block
pair given their total are Dirichlet-multinomial with concentration \eqn{k} per
edge; the Bayesian form puts a Beta(\code{alpha}, \code{beta}) prior on the
negative-binomial probability. At \eqn{k = 1} it equals
\code{"discrete-geometric"} exactly. Overdispersed counts fitted as
\code{"discrete-poisson"} invent block structure to explain the extra variance:
on random graphs without blocks and with counts of variance/mean 7-12, Poisson
fits found 8-31 blocks; negative binomial at the true \eqn{k} found 1. Choose
\eqn{k} by comparing description lengths over a grid of values (each value is a
separate model; add the log of the grid size for the choice).

All description lengths match graph-tool 3.9's \code{WeightedBlockState} and
\code{NestedBlockState(base_state=WeightedBlockState)} to machine precision
(including the exact continuous quantization scale
\eqn{\epsilon = 0.1 \min_{w_e \neq 0} |w_e|} and the nested
\code{real-normal} variance hyperprior).

Weighted states integrate directly with \code{\link{block_membership}},
\code{\link{block_sizes}}, \code{\link{nested_levels}},
\code{\link{nested_project}}, \code{\link{nested_summary}},
\code{\link{plot_hierarchy}}, \code{\link{mcmc_sweep}},
\code{\link{mcmc_equilibrate}}, \code{\link{mcmc_anneal}},
\code{\link{sample_posterior}}, and \code{\link{model_evidence}}.
}
\value{
Flat constructors and minimizers return a \code{rgraphtool_weighted_block_state};
nested constructors and minimizers return a
\code{rgraphtool_nested_weighted_block_state}. Description length functions return
a numeric scalar (or per-level numeric vector for
\code{nested_weighted_level_entropy()}).
}
\references{
Peixoto, T. P. (2018). Nonparametric weighted stochastic block models.
\emph{Physical Review E} 97, 012306.
}
\seealso{\code{\link{compare_weight_models}}, \code{\link{transform_weights}},
\code{\link{graph_set_edge_property}}, \code{demo("weighted_sbm")}}
\examples{
set.seed(42)
g <- graph_create(12, directed = FALSE)
edges <- rbind(t(utils::combn(1:6, 2)), t(utils::combn(7:12, 2)), c(6, 7))
w <- c(stats::rnorm(15, 5, 0.5), stats::rnorm(15, -5, 0.5), 0)
graph_add_edges(g, edges, properties = list(w = w))
st <- minimize_weighted_blockmodel_dl(g, weights = "w", types = "real-normal")
st
nst <- minimize_nested_weighted_blockmodel_dl(g, weights = "w",
                                              types = "real-normal")
nst
}
