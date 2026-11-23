\name{pp_block_state}
\alias{pp_block_state}
\alias{pp_description_length}
\alias{pp_assortativity}
\alias{minimize_pp_blockmodel_dl}
\alias{compare_assortative_models}
\alias{print.graphtool_pp_block_state}
\alias{print.graphtool_assortative_model_comparison}
\title{Assortative planted-partition stochastic block model}
\usage{
pp_block_state(graph, blocks = NULL, uniform = FALSE,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"))
pp_description_length(state, degree_prior = NULL)
pp_assortativity(state)
minimize_pp_blockmodel_dl(graph, uniform = FALSE, degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"), B_min = 1L, B_max = Inf,
  refine = TRUE, epsilon = 1e-3)
compare_assortative_models(graph,
  models = c("pp", "pp_uniform", "sbm_flat", "sbm_nested"),
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  method = c("map", "posterior"),
  evidence_method = c("variational", "importance"),
  nruns = 1L, nsamples = 200L, thin = 5L, burnin = NULL, mle = TRUE,
  B_min = 1L, B_max = Inf, refine = TRUE)
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
  \item{state}{A \code{graphtool_pp_block_state} created by
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
  \item{nsamples, thin, burnin, mle}{Posterior sampling and entropy estimation
  controls forwarded to \code{\link{model_evidence}}.}
}
\description{
Native C++ implementation of the microcanonical Bayesian planted-partition
stochastic block model (Zhang and Peixoto 2020), matching graph-tool's
\code{PPBlockState}. Constrains community structure to assortative (or
disassortative) block patterns parameterized only by within-community edges
\eqn{E_{\mathrm{in}}} and between-community edges \eqn{E_{\mathrm{out}}}, and
provides declarative model comparison against the general flat and nested SBMs.
}
\details{
Unlike the general SBM, which estimates an arbitrary \eqn{B \times B} block
affinity matrix, the planted-partition model restricts off-diagonal interactions
to a single global rate while modeling internal edges either per-community
(\code{uniform = FALSE}) or uniformly (\code{uniform = TRUE}). This gives
\eqn{O(B)} state memory and zero-allocation incremental moves while remaining
free of modularity's resolution limit and overfitting on random graphs.

Planted-partition states work directly with \code{\link{block_membership}},
\code{\link{block_sizes}}, \code{\link{mcmc_sweep}},
\code{\link{mcmc_equilibrate}}, \code{\link{mcmc_anneal}},
\code{\link{sample_posterior}}, and \code{\link{model_evidence}}.
}
\value{
\code{pp_block_state()} and \code{minimize_pp_blockmodel_dl()} return a
\code{graphtool_pp_block_state} list with \code{blocks}, \code{graph},
\code{uniform}, \code{degree_corrected}, \code{degree_prior}, \code{B}, and
optional \code{fit}.

\code{pp_description_length()} returns the exact microcanonical description
length in nats.

\code{pp_assortativity()} returns a list with \code{E_in}, \code{E_out},
\code{E_total}, \code{fraction_in}, \code{B}, effective group count \code{Be},
and per-group internal edge counts \code{err}.

\code{compare_assortative_models()} returns a
\code{graphtool_assortative_model_comparison} list with ranked \code{summary},
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
