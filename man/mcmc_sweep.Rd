\name{mcmc_sweep}
\alias{mcmc_sweep}
\alias{mcmc_equilibrate}
\alias{mcmc_anneal}
\alias{sample_posterior}
\title{MCMC posterior sampling, equilibration, and simulated annealing for SBMs}
\usage{
mcmc_sweep(state, beta = 1, niter = 1L, c = NULL, d = 0.01,
  multiflip = TRUE, B_min = 1L, B_max = Inf, psingle = NULL,
  psplit = 1, pmerge = 1, pmergesplit = 1, gibbs_sweeps = 10L,
  ls_shuffle = TRUE)
mcmc_equilibrate(state, wait = 1000L, nbreaks = 2L, max_niter = Inf,
  force_niter = NULL, epsilon = 0, multiflip = TRUE, mcmc_args = list(),
  history = FALSE, callback = NULL, collect = FALSE, thin = 1L,
  verbose = FALSE)
mcmc_anneal(state, beta_range = c(1, 10), niter = 100L,
  history = FALSE, mcmc_equilibrate_args = list(), verbose = FALSE)
sample_posterior(state, nsamples = 1000L, thin = 10L, burnin = NULL,
  multiflip = TRUE, mcmc_args = list(), equilibrate_args = list(),
  history = FALSE, verbose = FALSE)
}
\arguments{
  \item{state}{A flat block state from \code{\link{block_state}},
  \code{\link{minimize_blockmodel_dl}}, \code{\link{weighted_block_state}},
  \code{\link{minimize_weighted_blockmodel_dl}},
  \code{\link{layered_block_state}},
  \code{\link{minimize_layered_blockmodel_dl}},
  \code{\link{pp_block_state}}, \code{\link{minimize_pp_blockmodel_dl}},
  \code{\link{ranked_block_state}}, or
  \code{\link{minimize_ranked_blockmodel_dl}}, or a nested block state from
  \code{\link{nested_block_state}},
  \code{\link{minimize_nested_blockmodel_dl}},
  \code{\link{nested_weighted_block_state}},
  \code{\link{minimize_nested_weighted_blockmodel_dl}},
  \code{\link{nested_layered_block_state}},
  \code{\link{minimize_nested_layered_blockmodel_dl}},
  \code{\link{nested_ranked_block_state}}, or
  \code{\link{minimize_nested_ranked_blockmodel_dl}}.}
  \item{beta}{Inverse temperature (\eqn{\beta \ge 0}, or \code{Inf} for greedy
  zero-temperature moves). Posterior sampling uses \code{beta = 1}.}
  \item{niter}{Number of sweeps per call (for \code{mcmc_sweep}), or number of
  exponential cooling steps (for \code{mcmc_anneal}).}
  \item{c}{Ergodic proposal smoothing parameter (\eqn{c \ge 0}; \code{Inf} gives
  uniform proposals). Defaults to \code{0.5} for flat states and \code{1}
  (scaled by \eqn{2^l} at level \eqn{l}) for nested states. Nested states also
  accept a numeric vector of per-level values.}
  \item{d}{Probability in \eqn{[0, 1)} of proposing a new (empty) group during a
  single-node move.}
  \item{multiflip}{If \code{TRUE} (the default), use merge-split MCMC combining
  single-node moves, group splits, group merges, and merge-split moves; if
  \code{FALSE}, use single-node Metropolis-Hastings sweeps only.}
  \item{B_min, B_max}{Bounds on the number of occupied groups (applied to level 0
  for nested states).}
  \item{psingle, psplit, pmerge, pmergesplit}{Relative move-type weights for
  merge-split MCMC. \code{psingle = NULL} defaults to the number of nodes at the
  current level.}
  \item{gibbs_sweeps}{Number of Gibbs sweeps used to refine split and merge-split
  proposals (at least 1).}
  \item{ls_shuffle}{For nested states, whether to randomize the order in which
  hierarchy levels are swept.}
  \item{wait, nbreaks, max_niter, force_niter, epsilon}{Equilibration controls
  matching graph-tool's \code{mcmc_equilibrate()}: waits for \code{nbreaks}
  consecutive intervals of \code{wait} iterations without a record-breaking
  change in description length, up to \code{max_niter} iterations (or runs
  exactly \code{force_niter} iterations when specified).}
  \item{mcmc_args}{Named list of arguments forwarded to \code{mcmc_sweep()}.}
  \item{history}{Whether to return a data frame recording \code{nattempts},
  \code{nmoves}, \code{entropy}, and \code{beta} at each iteration.}
  \item{callback}{Optional R function called as \code{callback(state)} after each
  iteration. When \code{NULL}, the entire equilibration loop runs in C++.}
  \item{collect}{Whether to collect sampled partitions (every \code{thin}
  iterations) during equilibration.}
  \item{thin}{Thinning interval (number of equilibration steps between collected
  samples).}
  \item{verbose}{Whether to print iteration progress.}
  \item{beta_range}{Numeric vector of length 2 giving the starting and ending
  inverse temperatures for exponential cooling.}
  \item{mcmc_equilibrate_args}{Named list of equilibration arguments used at each
  annealing step.}
  \item{nsamples}{Number of thinned posterior partitions to collect.}
  \item{burnin}{Number of burn-in sweeps before sampling. When \code{NULL} (the
  default), runs \code{mcmc_equilibrate()} with \code{equilibrate_args} first.}
  \item{equilibrate_args}{Named list of equilibration controls used when
  \code{burnin = NULL}.}
}
\description{
Native C++ Metropolis-Hastings and merge-split Markov chain Monte Carlo samplers
for flat and nested microcanonical stochastic block models, targeting the
exact posterior \eqn{P(\mathbf{b} \mid A) \propto \exp(-\beta \Sigma(\mathbf{b}))}.
These are the R counterparts of graph-tool's \code{mcmc_sweep()},
\code{multiflip_mcmc_sweep()}, \code{mcmc_equilibrate()}, and
\code{mcmc_anneal()}.
}
\details{
Single-node moves use Peixoto's neighbor-based ergodic proposal with exact
forward and reverse Hastings probabilities. Merge-split sweeps combine
single-node moves with group splits, group merges, and merge-split moves refined
by \code{gibbs_sweeps} conditional Gibbs sweeps using an auxiliary-variable
Hastings ratio that preserves exact detailed balance for any
\code{gibbs_sweeps >= 1} and any finite \eqn{\beta \ge 0}.

For nested states at level \eqn{l < L - 1}, new groups created by single-node
moves or splits inherit the parent label at level \eqn{l + 1}; correspondingly,
emptying a singleton group or merging two groups at level \eqn{l} is restricted
to groups sharing the same parent, while moves at level \eqn{l + 1} reassign
groups across parents. This maintains ergodicity and exact detailed balance over
unlabelled hierarchies.

Unlike Python's mutating methods, these R functions treat states as immutable and
return the updated state in the \code{$state} element of the result. All random
draws use R's RNG, so \code{set.seed()} reproduces trajectories. As with any
MCMC sampler on multimodal posteriors, record-breaking equilibration is a
heuristic and does not guarantee global mixing across well-separated modes.
}
\value{
\code{mcmc_sweep()}, \code{mcmc_equilibrate()}, and \code{mcmc_anneal()} return
a list with the updated \code{state}, final \code{entropy}, cumulative \code{dS},
\code{nattempts}, \code{nmoves}, \code{niter}, and optionally \code{history} and
\code{samples}. \code{sample_posterior()} returns a list with final \code{state},
\code{samples} (list of 1-based partitions or nested hierarchies), \code{nested},
\code{equilibration}, and \code{sampling}.
}
\references{
Peixoto, T. P. (2014). Efficient Monte Carlo and greedy heuristic for the
inference of stochastic block models. \emph{Physical Review E} 89, 012804.
doi:10.1103/PhysRevE.89.012804

Peixoto, T. P. (2020). Merge-split Markov chain Monte Carlo for community
detection. \emph{Physical Review E} 102, 012305.
doi:10.1103/PhysRevE.102.012305
}
\seealso{\code{\link{partition_mode_state}}, \code{\link{mode_cluster_state}},
\code{\link{plot_marginals}}, \code{demo("posterior_sampling")}}
\examples{
g <- graph_create(12, directed = FALSE)
graph_add_edges(g, rbind(t(combn(1:6, 2)), t(combn(7:12, 2)), c(1, 7)))
set.seed(1)
fit <- minimize_blockmodel_dl(g)
post <- sample_posterior(fit, nsamples = 50, thin = 2, burnin = 20)
pm <- partition_mode_state(post)
pm
}
