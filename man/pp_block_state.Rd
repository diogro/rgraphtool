\name{pp_block_state}
\alias{pp_block_state}
\alias{pp_description_length}
\alias{pp_assortativity}
\alias{minimize_pp_blockmodel_dl}
\alias{compare_assortative_models}
\alias{print.rgraphtool_pp_block_state}
\alias{print.rgraphtool_assortative_model_comparison}
\title{Assortative planted-partition stochastic block model}
\usage{
pp_block_state(graph, blocks = NULL, uniform = FALSE,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  weights = NULL, types = NULL, params = "microcanonical",
  layers = NULL, independent = FALSE)
pp_description_length(state, degree_prior = NULL)
pp_assortativity(state)
minimize_pp_blockmodel_dl(graph, uniform = FALSE, degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"), B_min = 1L, B_max = Inf,
  refine = TRUE, epsilon = 1e-3,
  weights = NULL, types = NULL, params = "microcanonical",
  layers = NULL, independent = FALSE)
compare_assortative_models(graph,
  models = c("pp", "pp_uniform", "sbm_flat", "sbm_nested"),
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  method = c("map", "posterior"),
  evidence_method = c("variational", "importance"),
  nruns = 1L, cores = NULL, nsamples = 200L, thin = 5L, burnin = NULL, mle = TRUE,
  B_min = 1L, B_max = Inf, refine = TRUE,
  weights = NULL, types = NULL, params = "microcanonical")
}
\arguments{
  \item{graph}{A graph created by \code{\link{graph_create}}.}
  \item{blocks}{Optional integer vector of 1-based group memberships of length
  \code{graph_num_vertices(graph)}. Defaults to a single group.}
  \item{uniform}{If \code{FALSE} (the default), use the non-uniform
  planted-partition model where each community \eqn{r} has its own internal edge
  count \eqn{e_{rr}}; if \code{TRUE}, constrain all communities to share a
  uniform internal rate.}
  \item{degree_corrected}{Whether to fit the degree-corrected (\code{TRUE}) or
  non-degree-corrected (\code{FALSE}) planted-partition model.}
  \item{degree_prior}{Degree sequence hyperprior (\code{"distributed"} or
  \code{"uniform"}).}
  \item{state}{A \code{rgraphtool_pp_block_state} created by
  \code{pp_block_state()} or \code{minimize_pp_blockmodel_dl()}.}
  \item{B_min, B_max}{Minimum and maximum number of occupied groups.}
  \item{refine}{Whether to run greedy refinement sweeps after merge-split
  minimization.}
  \item{epsilon}{Stopping tolerance in nats for greedy refinement.}
  \item{models}{Character vector naming a subset of \code{"pp"},
  \code{"pp_uniform"}, \code{"sbm_flat"}, and \code{"sbm_nested"} to compare.}
  \item{method}{Compare models by \code{"map"} minimum description length or
  full Bayesian \code{"posterior"} marginal evidence.}
  \item{evidence_method}{Evidence estimator forwarded to
  \code{\link{model_evidence}} when \code{method = "posterior"}.}
  \item{nruns}{Number of independent minimization restarts per model.}
  \item{cores}{Number of worker processes running the \code{nruns} independent starts in parallel (forked; \code{NULL} uses \code{getOption("rgraphtool.cores", 1L)}). Workers run single-threaded; results do not depend on \code{cores}. Unavailable on Windows.}
  \item{nsamples, thin, burnin, mle}{Posterior sampling and entropy estimation
  controls forwarded to \code{\link{model_evidence}}.}
  \item{weights}{Optional edge covariate specification (edge property name(s),
  vector of length \code{graph_num_edges(graph)}, or matrix/data.frame with one
  column per covariate) using the same distributions and conventions as
  \code{\link{weighted_block_state}}. When provided, \code{compare_assortative_models()}
  compares weighted planted-partition models against weighted flat and nested
  SBMs (\code{\link{minimize_weighted_blockmodel_dl}} and
  \code{\link{minimize_nested_weighted_blockmodel_dl}}).}
  \item{types}{Optional character vector of covariate families corresponding to
  \code{weights} (\code{"real-normal"}, \code{"real-exponential"},
  \code{"discrete-poisson"}, \code{"discrete-geometric"}, or
  \code{"discrete-binomial"}). Required when \code{weights} is supplied.}
  \item{params}{Hyperprior specification for \code{weights}: \code{"microcanonical"}
  (default) or explicit Bayesian prior hyperparameters as in
  \code{\link{weighted_block_state}}.}
  \item{layers}{Optional discrete edge layer labels: an edge property name or a
  vector of length \code{graph_num_edges(graph)} (as in
  \code{\link{layered_block_state}}). Enables flat multi-layer planted-partition
  inference, optionally combined with \code{weights}.}
  \item{independent}{Used when \code{layers} is provided: if \code{FALSE}
  (default), fit the discrete edge-covariate PP model with pooled degrees and
  pooled covariate buckets across layers; if \code{TRUE}, fit conditionally
  independent PP layers with layer-specific degrees, active group sizes, and
  (when weighted) per-layer covariate distributions.}
}
\description{
Native C++ implementation of the microcanonical Bayesian planted-partition
stochastic block model (Zhang and Peixoto 2020), matching graph-tool's
\code{PPBlockState}. Constrains community structure to assortative (or
disassortative) block patterns parameterized only by within-community edges
\eqn{E_{\mathrm{in}}} and between-community edges \eqn{E_{\mathrm{out}}},
optionally paired with exponential-family or microcanonical edge-weight
covariate distributions, and provides declarative model comparison against the
general flat and nested SBMs.
}
\details{
Unlike the general SBM, which estimates an arbitrary \eqn{B \times B} block
affinity matrix, the planted-partition model restricts off-diagonal interactions
to a single global rate while modeling internal edges either per-community
(\code{uniform = FALSE}) or uniformly (\code{uniform = TRUE}). This gives
\eqn{O(B)} state memory and zero-allocation incremental moves while remaining
free of modularity's resolution limit and overfitting on random graphs.

\strong{Weighted planted-partition models and caveats:}
When \code{weights} and \code{types} are supplied, each active edge bucket
receives its own weight distribution under the exact conjugate/microcanonical
likelihoods of \code{\link{weighted_block_state}} (2 buckets \code{in}/\code{out}
when \code{uniform = TRUE}; \eqn{B + 1} buckets \eqn{\{rr\}} and \code{out} when
\code{uniform = FALSE}):
\itemize{
  \item \emph{Pooled inter-module background:} All between-community edge
  weights are pooled into a single \code{out} distribution rather than
  modeling pair-specific \eqn{(r, s)} block weight differences. If specific
  pairs of modules have distinct positive vs. negative weight patterns (e.g.,
  inter-module antagonism in a gene regulatory network), use
  \code{compare_assortative_models(..., weights = ..., types = ...)} to test
  whether \code{"sbm_flat"} or \code{"sbm_nested"} has higher Bayesian evidence.
  \item \emph{Mean-ordering / sign symmetry (especially at \eqn{B = 2}):}
  With signed covariates such as \code{"real-normal"}, the covariate likelihood
  tests whether the within-community and between-community weight distributions
  \emph{differ}, without enforcing \eqn{\mu_{\mathrm{in}} > \mu_{\mathrm{out}}}.
  On strongly anti-correlated bipartitions (particularly at \eqn{B = 2}), inspect
  \code{pp_assortativity()} (\code{fraction_in}, \code{weight_in_mean},
  \code{weight_out_mean}) to confirm whether a recovered partition is assortative
  or disassortative.
  \item \emph{Extension beyond upstream graph-tool 3.9:} Upstream
  \code{graph-tool} 3.9 \code{PPBlockState} does not support edge covariates
  (\code{recs}) or edge layers; \code{rgraphtool} extends PP to both
  (flat, single-level) layers and covariates using the exact \code{graph-tool}
  microcanonical/conjugate likelihoods. At \eqn{B = 1}, the weighted/layered
  PP description length matches \code{\link{weighted_block_state}} or
  \code{\link{layered_block_state}} identically, and \code{compare_layered_models()}
  supports direct comparison against general layered SBMs.
}

Planted-partition states work directly with \code{\link{block_membership}},
\code{\link{block_sizes}}, \code{\link{mcmc_sweep}},
\code{\link{mcmc_equilibrate}}, \code{\link{mcmc_anneal}},
\code{\link{sample_posterior}}, and \code{\link{model_evidence}}.
}
\value{
\code{pp_block_state()} and \code{minimize_pp_blockmodel_dl()} return a
\code{rgraphtool_pp_block_state} list with \code{blocks}, \code{graph},
\code{uniform}, \code{degree_corrected}, \code{degree_prior}, \code{B}, and
optional \code{fit} (plus \code{weights}, \code{types}, \code{params},
\code{weighted}, \code{layers}, \code{layer_levels}, \code{C},
\code{independent}, and \code{layered}).

\code{pp_description_length()} returns the exact microcanonical description
length in nats.

\code{pp_assortativity()} returns a list with \code{E_in}, \code{E_out},
\code{E_total}, \code{fraction_in}, \code{B}, effective group count \code{Be},
and per-group internal edge counts \code{err} (plus \code{weight_in_mean} and
\code{weight_out_mean} when \code{state} is weighted, and per-layer fields
\code{E_in_layer}, \code{E_out_layer}, \code{E_total_layer},
\code{fraction_in_layer}, \code{err_layer}, \code{weight_in_mean_layer}, and
\code{weight_out_mean_layer} when \code{state} is layered).

\code{compare_assortative_models()} returns a
\code{rgraphtool_assortative_model_comparison} list with ranked \code{summary},
\code{best_model}, \code{best_state}, \code{states}, optional \code{evidences},
and \code{log_bayes_factors}.
}
\references{
Zhang, L. and Peixoto, T. P. (2020). Statistical inference of assortative
community structures. \emph{Physical Review Research} 2, 043271.
doi:10.1103/PhysRevResearch.2.043271
}
\seealso{\code{\link{minimize_blockmodel_dl}},
\code{\link{minimize_nested_blockmodel_dl}}, \code{\link{mcmc_sweep}},
\code{\link{model_evidence}}, \code{demo("assortative_sbm")}}
\examples{
g <- graph_create(12, directed = FALSE)
graph_add_edges(g, rbind(t(combn(1:6, 2)), t(combn(7:12, 2)), c(1, 7)))
set.seed(1)
st <- minimize_pp_blockmodel_dl(g)
st
pp_assortativity(st)
}
