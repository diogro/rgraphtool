\name{compare_layered_models}
\alias{compare_layered_models}
\alias{print.rgraphtool_layered_model_comparison}
\title{Compare edge-covariate and independent-layers SBM variants}
\usage{
compare_layered_models(graph, layers,
  models = c("edge_covariates", "independent_layers"), nested = TRUE,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  method = c("map", "posterior"),
  evidence_method = c("variational", "importance"), nruns = 1L,
  nsamples = 200L, thin = 5L, burnin = NULL, mle = TRUE, B_min = 1L,
  B_max = Inf, refine = TRUE)
}
\arguments{
  \item{graph}{A graph created by \code{\link{graph_create}}.}
  \item{layers}{Discrete edge layer labels: an edge property name or a vector
  of length \code{graph_num_edges(graph)}.}
  \item{models}{Character vector naming a subset of \code{"edge_covariates"},
  \code{"independent_layers"}, \code{"edge_covariates_dc"},
  \code{"edge_covariates_ndc"}, \code{"independent_layers_dc"}, and
  \code{"independent_layers_ndc"}. The unqualified names use
  \code{degree_corrected}.}
  \item{nested}{Whether to fit nested (\code{TRUE}, default) or flat
  (\code{FALSE}) layered SBMs.}
  \item{degree_corrected}{Degree correction setting used for
  \code{"edge_covariates"} and \code{"independent_layers"}.}
  \item{degree_prior}{Degree sequence hyperprior (\code{"distributed"} or
  \code{"uniform"}).}
  \item{method}{Rank candidate models by \code{"map"} minimum description length
  or full Bayesian \code{"posterior"} marginal evidence.}
  \item{evidence_method}{Evidence estimator forwarded to
  \code{\link{model_evidence}} when \code{method = "posterior"}.}
  \item{nruns}{Number of independent minimization restarts per candidate.}
  \item{nsamples, thin, burnin, mle}{Posterior sampling arguments used when
  \code{method = "posterior"}.}
  \item{B_min, B_max}{Bounds on the number of level-0 groups.}
  \item{refine}{Whether to run local refinement sweeps during minimization.}
}
\description{
Declarative, callback-free model selection between the discrete edge-covariate
layered SBM (\code{independent = FALSE}) and the conditionally independent
layers SBM (\code{independent = TRUE}), optionally across degree-corrected and
non-degree-corrected variants and using either MAP description lengths or full
Bayesian posterior marginal evidence.
}
\value{
A \code{rgraphtool_layered_model_comparison} list containing:
\itemize{
  \item \code{summary}: ranked data frame with \code{model}, \code{independent},
  \code{degree_corrected}, \code{B}, \code{description_length},
  \code{posterior_entropy}, \code{neg_log_evidence}, \code{delta}, and
  \code{posterior_weight}.
  \item \code{best_model}: name of the winning layered model.
  \item \code{best_state}: fitted layered state for the winning model.
  \item \code{states}: named list of fitted states ordered from best to worst.
  \item \code{evidences}: named list of \code{rgraphtool_model_evidence} objects
  (when \code{method = "posterior"}).
  \item \code{log_bayes_factors}: pairwise matrix of log Bayes factors.
}
}
\references{
Peixoto, T. P. (2015). Inferring the mesoscale structure of layered,
edge-valued, and time-varying networks. \emph{Physical Review E} 92, 042807.
}
\seealso{\code{\link{layered_block_state}}, \code{\link{model_evidence}},
\code{demo("layered_sbm")}}
\examples{
set.seed(42)
g <- graph_create(12, directed = FALSE)
edges <- rbind(t(utils::combn(1:6, 2)), t(utils::combn(7:12, 2)))
lay <- c(rep("layer_a", 15), rep("layer_b", 15))
graph_add_edges(g, edges)
cmp <- compare_layered_models(g, layers = lay, nested = FALSE)
cmp
}
