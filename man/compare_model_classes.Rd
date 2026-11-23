\name{compare_model_classes}
\alias{compare_model_classes}
\alias{model_evidence}
\alias{print.graphtool_model_evidence}
\alias{print.graphtool_model_class_comparison}
\title{Posterior model evidence and SBM model class selection}
\usage{
model_evidence(state, nsamples = 200L, thin = 5L, burnin = NULL,
  beta = 1, mle = TRUE, method = c("variational", "importance"),
  samples = NULL, entropies = NULL, multiflip = TRUE,
  c = 0.5, d = 0.01, B_min = 1L, B_max = Inf)
compare_model_classes(graph,
  models = c("flat_ndc", "flat_dc", "nested_ndc", "nested_dc",
    "overlap_flat_ndc", "overlap_flat_dc",
    "overlap_nested_ndc", "overlap_nested_dc"),
  degree_prior = c("distributed", "uniform"),
  method = c("map", "posterior"),
  evidence_method = c("variational", "importance"),
  nruns = 1L, nsamples = 200L, thin = 5L, burnin = NULL, mle = TRUE,
  B_min = 1L, B_max = Inf, refine = TRUE)
}
\arguments{
  \item{state}{A fitted \code{graphtool_block_state},
  \code{graphtool_nested_block_state}, \code{graphtool_overlap_block_state},
  \code{graphtool_nested_overlap_block_state},
  \code{graphtool_weighted_block_state},
  \code{graphtool_nested_weighted_block_state},
  \code{graphtool_layered_block_state},
  \code{graphtool_nested_layered_block_state},
  \code{graphtool_pp_block_state},
  \code{graphtool_ranked_block_state}, or
  \code{graphtool_nested_ranked_block_state} object.}
  \item{graph}{A graph created by \code{\link{graph_create}}.}
  \item{models}{Character vector naming a subset of the 8 SBM model classes:
  \code{"flat_ndc"}, \code{"flat_dc"}, \code{"nested_ndc"}, \code{"nested_dc"},
  \code{"overlap_flat_ndc"}, \code{"overlap_flat_dc"},
  \code{"overlap_nested_ndc"}, and \code{"overlap_nested_dc"}.}
  \item{degree_prior}{Degree sequence hyperprior (\code{"distributed"} or
  \code{"uniform"}).}
  \item{method}{For \code{compare_model_classes()}, rank models by \code{"map"}
  minimum description length or full Bayesian \code{"posterior"} marginal
  evidence. For \code{model_evidence()}, choose between \code{"variational"}
  (\eqn{\ln P(A) = -\langle \Sigma \rangle + H(b)}) and \code{"importance"}
  (Peixoto 2021 Eq. 29 mode-mixture importance sampling).}
  \item{evidence_method}{Evidence estimator forwarded to \code{model_evidence()}
  when \code{method = "posterior"}.}
  \item{nruns}{Number of independent minimization restarts per model class,
  retaining the lowest description length fit.}
  \item{nsamples}{Number of thinned posterior samples to collect.}
  \item{thin}{Number of MCMC sweeps between collected samples.}
  \item{burnin}{Number of burn-in sweeps before sampling, or \code{NULL} to
  equilibrate automatically (flat/nested) or use 30 sweeps (overlap).}
  \item{beta}{Inverse temperature (default \code{1}).}
  \item{mle}{Whether to use MLE or Dirichlet(1) posterior mean marginals for
  posterior entropy and mode log-probabilities.}
  \item{samples, entropies}{Optional pre-collected posterior partition samples
  and matching description lengths (e.g. from \code{\link{sample_posterior}}),
  avoiding redundant MCMC runs.}
  \item{multiflip}{Whether to use merge-split MCMC moves for flat/nested states.}
  \item{c, d}{MCMC proposal mixing and empty-group parameters.}
  \item{B_min, B_max}{Minimum and maximum number of level-0 groups.}
  \item{refine}{Whether to run refinement sweeps during minimization.}
}
\description{
Declarative, callback-free Bayesian model class selection across flat,
hierarchical (nested), degree-corrected, non-degree-corrected, and overlapping
stochastic block models. Computes both MAP description lengths and full
posterior marginal log-evidences \eqn{\ln P(A | \mathcal{M})}, posterior model
weights, and pairwise log Bayes factors.
}
\details{
For a given model class, \code{model_evidence()} equilibrates an MCMC chain,
collects posterior partition samples and exact description lengths
\eqn{\Sigma(b^{(m)})} natively in C++, clusters the samples into posterior modes
via \code{\link{mode_cluster_state}}, and computes both the variational
evidence \eqn{\ln P(A) = -\langle \Sigma \rangle + H_{\mathrm{posterior}}}
(matching graph-tool's model class selection cookbook) and Peixoto's (2021)
mode-mixture importance estimate
\eqn{-\ln P(A) = \ln\!\left(\frac{1}{M}\sum_{m=1}^M \exp(\Sigma(b^{(m)}) + \ln q(b^{(m)}))\right)}.

Note that for overlapping states, the posterior partition entropy is estimated
from per-vertex majority-group samples rather than full half-edge trajectories,
providing a conservative lower bound on the half-edge posterior entropy.
}
\value{
\code{model_evidence()} returns a \code{graphtool_model_evidence} list with
\code{log_evidence}, \code{neg_log_evidence}, \code{variational_log_evidence},
\code{importance_log_evidence}, \code{mean_dl}, \code{sd_dl}, \code{map_dl},
\code{posterior_entropy}, \code{cluster}, \code{samples}, \code{entropies}, and
\code{state}.

\code{compare_model_classes()} returns a \code{graphtool_model_class_comparison}
list with ranked \code{summary} data frame, \code{best_model},
\code{best_state}, named list of fitted \code{states}, optional
\code{evidences}, and pairwise \code{log_bayes_factors} matrix.
}
\references{
Peixoto, T. P. (2015). Model selection and hypothesis testing for large-scale
network models with overlapping groups. \emph{Physical Review X} 5, 011033.

Peixoto, T. P. (2021). Revealing consensus and dissensus between network
partitions. \emph{Physical Review X} 11, 021003.
}
\seealso{\code{\link{overlap_block_state}}, \code{\link{mode_cluster_state}},
\code{\link{sample_posterior}}, \code{demo("model_class_selection")}}
\examples{
g <- graph_create(n = 12L, directed = FALSE)
edges <- rbind(
  t(utils::combn(1:6, 2)),
  t(utils::combn(7:12, 2)),
  c(6, 7)
)
graph_add_edges(g, edges)
set.seed(42)
cmp <- compare_model_classes(g, models = c("flat_ndc", "flat_dc",
                                           "nested_ndc", "overlap_flat_ndc"))
cmp
}
