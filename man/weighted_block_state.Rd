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
  refine = TRUE, epsilon = 1e-3, nthreads = NULL)
nested_weighted_block_state(graph, weights, types, params = "microcanonical",
  levels = NULL, degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"))
nested_weighted_level_entropy(state, degree_prior = NULL)
nested_weighted_description_length(state, degree_prior = NULL)
minimize_nested_weighted_blockmodel_dl(graph, weights, types,
  params = "microcanonical", degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"), B_min = 1L, B_max = Inf,
  refine = TRUE, epsilon = 1e-3, nthreads = NULL)
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
  or \code{"discrete-binomial"}.}
  \item{params}{Covariate hyperprior specification: \code{"microcanonical"}
  (default, Peixoto 2018 maximum-entropy microcanonical priors),
  \code{"bayesian"} (empirical-Bayes conjugate hyperpriors matching graph-tool),
  a named list of hyperparameters (e.g. \code{list(alpha = 1, beta = 2)} or
  \code{list(N = 10, alpha = 1, beta = 1)} or
  \code{list(m0 = 0, k0 = 1, v0 = 1, nu0 = 3)}), or a list of length \eqn{K}
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
