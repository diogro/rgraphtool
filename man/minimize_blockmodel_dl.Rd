\name{minimize_blockmodel_dl}
\alias{minimize_blockmodel_dl}
\alias{refine_blockmodel_dl}
\alias{compare_blockmodels}
\alias{print.rgraphtool_blockmodel_comparison}
\title{Fit a flat stochastic block model by minimum description length}
\usage{
minimize_blockmodel_dl(graph, degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"), B_min = 1, B_max = Inf,
  refine = TRUE, epsilon = 1e-3, init_blocks = NULL, nthreads = NULL,
  control = list())
refine_blockmodel_dl(graph, state, degree_prior = NULL, epsilon = 1e-3,
  nthreads = NULL, control = list())
compare_blockmodels(graph, degree_prior = c("distributed", "uniform"),
  nested = FALSE, nruns = 1L, cores = NULL, ...)
}
\arguments{
  \item{graph}{A graph created by \code{graph_create()}.}
  \item{degree_corrected}{Whether to fit the degree-corrected model.}
  \item{degree_prior}{Degree-sequence prior, \code{"distributed"} (graph-tool's
  default) or \code{"uniform"}. For \code{refine_blockmodel_dl()}, \code{NULL}
  uses the prior stored in \code{state}.}
  \item{B_min, B_max}{Bounds on the number of nonempty groups.}
  \item{refine}{Whether to run the refinement loop after the multilevel search.}
  \item{epsilon}{Refinement stops when a pass changes the description length by a
  relative amount of at most \code{epsilon}.}
  \item{init_blocks}{Optional warm-start partition (either a 1-based integer
  vector of length \code{graph_num_vertices(graph)} or an existing
  \code{rgraphtool_block_state}), skipping initialization from singleton groups.}
  \item{nthreads}{Optional positive integer number of parallel worker threads
  (or \code{"auto"}). \code{NULL} uses \code{\link{rgraphtool_num_threads}}.}
  \item{control}{Named list of search settings; see Details.}
  \item{state}{A state created by \code{block_state()} or a fitting function.}
  \item{nested}{Whether \code{compare_blockmodels()} fits nested models with
  \code{\link{minimize_nested_blockmodel_dl}} instead of flat ones.}
  \item{nruns}{Number of fits per model; the one with the smallest description
  length is kept.}
  \item{cores}{Number of worker processes running the \code{nruns} independent starts in parallel (forked; \code{NULL} uses \code{getOption("rgraphtool.cores", 1L)}). Workers run single-threaded; results do not depend on \code{cores}. Unavailable on Windows.}
  \item{...}{Passed to the fitting function.}
}
\description{
\code{minimize_blockmodel_dl()} is the R counterpart of graph-tool's
\code{minimize_blockmodel_dl()} for the flat \code{BlockState}. It runs a native
multilevel agglomerative search followed by refinement, and returns the
partition with the smallest description length it found.
}
\details{
The search follows graph-tool 3.9's \code{multilevel_mcmc_sweep()} at zero
temperature (\eqn{\beta = \infty}) over all groups:
\enumerate{
  \item Every vertex starts in its own group. Single-vertex sweeps then move
  vertices into neighbouring groups, starting at \eqn{\beta =} \code{init_beta}
  for \code{init_min_iter_beta} sweeps and continuing at zero temperature for at
  least \code{init_min_iter} sweeps, until a sweep reduces the number of groups
  by less than a factor of \code{init_r}.
  \item A bisection search with random midpoints looks for the number of groups
  \eqn{B} that minimizes the description length. To reach a smaller \eqn{B}
  from a cached partition, the search alternates merges (each group tries
  \code{merge_sweeps} neighbour-proposed merge partners, and the best merges are
  applied until \eqn{B} falls by a factor of \code{r}) with up to
  \code{mh_sweeps} zero-temperature single-vertex sweeps. Every partition visited
  is cached by \eqn{B}.
  \item Refinement repeatedly merges groups when that lowers the description
  length and sweeps vertices among the existing groups. It stops after a pass
  whose relative change is at most \code{epsilon}.
}
Moves never increase the description length after the initial
\code{init_beta} sweeps. The result is a local optimum, not a guaranteed
global one, so run several seeds and keep the smallest description length, as
the graph-tool documentation recommends.

\code{control} accepts \code{r} (0.9), \code{merge_sweeps} (10),
\code{mh_sweeps} (10), \code{init_r} (0.99), \code{init_min_iter} (8),
\code{init_min_iter_beta} (4), \code{init_beta} (1), \code{c} (0.5, the random
proposal weight used when a group has no neighbour-proposed merge partner),
\code{max_refine} (1000), \code{bisection} (\code{TRUE}; set to \code{FALSE} for
strictly agglomerative geometric coarsening without bisection backtracking on
large networks), \code{nthreads} (overrides \code{\link{rgraphtool_num_threads}}
when \code{nthreads = NULL}), and \code{verbose} (\code{FALSE}).

When \code{nthreads > 1}, candidate merge evaluations in \code{merge_sweep()}
and zero-temperature vertex proposal screening in \code{mh_sweep()} run in
parallel across a portable C++23 thread pool, while all state mutations and
priority-queue merge commits execute deterministically on the main thread. For
very large networks (\eqn{N \ge 10^4}), combining \code{nthreads} with strictly
agglomerative search (\code{control = list(bisection = FALSE, r = 0.5)}) keeps
only a single snapshot in memory and reduces runtime by an order of magnitude.

\code{compare_blockmodels()} fits the non-degree-corrected (NDC) and
degree-corrected (DC) models and reports which has the smaller description
length \eqn{\Sigma}. Both description lengths include all priors, so with
equal prior probabilities for the two models the posterior odds ratio is
\deqn{\Lambda = \frac{P(b, H_{NDC} \mid A)}{P(b, H_{DC} \mid A)}
  = e^{-(\Sigma_{NDC} - \Sigma_{DC})},}
as in the graph-tool cookbook. \code{log_posterior_odds} is
\eqn{\ln \Lambda = \Sigma_{DC} - \Sigma_{NDC}}, so negative values favour the
degree-corrected model. (The cookbook prints this difference with the opposite
sign for the football network.) Because the fits are stochastic local optima,
increase \code{nruns} before relying on small differences.
}
\value{
\code{minimize_blockmodel_dl()} and \code{refine_blockmodel_dl()} return a
block state (see \code{\link{block_state}}) with labels \code{1, ..., B}. Its
\code{fit} element holds \code{description_length}. For
\code{minimize_blockmodel_dl()} it also holds
\code{multilevel_description_length}, the score before refinement, and
\code{trace}, a data frame of every cached \eqn{(B, \mathrm{DL})} pair in the
order visited. The description length of \code{refine_blockmodel_dl()}'s result
is never higher than that of the input state.

\code{compare_blockmodels()} returns a list with \code{comparison}, a data frame
of models and description lengths; \code{log_posterior_odds};
\code{best_model}; \code{best_state}; \code{states}, both fitted states; and
\code{nested}. Printing it shows the description lengths and \eqn{\ln \Lambda}.
}
\references{
Peixoto, T. P. (2014). Efficient Monte Carlo and greedy heuristic for the
inference of stochastic block models. \emph{Physical Review E} 89, 012804.
doi:10.1103/PhysRevE.89.012804

Peixoto, T. P. (2017). Nonparametric Bayesian inference of the microcanonical
stochastic block model. \emph{Physical Review E} 95, 012317.
doi:10.1103/PhysRevE.95.012317
}
\seealso{\code{\link{block_state}}}
\examples{
g <- graph_create(12, directed = FALSE)
graph_add_edges(g, do.call(rbind, lapply(list(1:6, 7:12),
                                         function(v) t(combn(v, 2)))))
set.seed(1)
fit <- minimize_blockmodel_dl(g)
block_membership(fit)
block_description_length(fit)
head(fit$fit$trace)
compare_blockmodels(g)$comparison
}
