\name{uncertain_block_state}
\alias{measured_block_state}
\alias{mixed_measured_block_state}
\alias{uncertain_block_state}
\alias{latent_multigraph_block_state}
\alias{uncertain_description_length}
\alias{measured_error_posterior}
\alias{edge_probabilities}
\alias{reconstruct_network}
\alias{latent_multigraph}
\alias{print.rgraphtool_uncertain_block_state}
\alias{print.rgraphtool_reconstruction}
\alias{print.rgraphtool_latent_multigraph}
\title{Uncertain network reconstruction and latent Poisson multigraphs}
\usage{
measured_block_state(graph, n = NULL, x = NULL,
  n_default = 1L, x_default = 0L,
  fn_params = list(alpha = 1, beta = 1),
  fp_params = list(mu = 1, nu = 1),
  lp = NA_real_, lq = NA_real_,
  nested = TRUE, blocks = NULL,
  latent_edges = NULL, latent_weights = NULL,
  degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"),
  self_loops = FALSE, q_max = Inf)

mixed_measured_block_state(graph, n = NULL, x = NULL,
  n_default = 1L, x_default = 0L,
  fn_params = list(alpha = 1, beta = 1),
  fp_params = list(mu = 1, nu = 1),
  nested = TRUE, blocks = NULL,
  latent_edges = NULL, latent_weights = NULL,
  degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"),
  self_loops = FALSE, q_max = Inf)

uncertain_block_state(graph, q = NULL, q_default = 0,
  nested = TRUE, blocks = NULL,
  latent_edges = NULL, latent_weights = NULL,
  degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"),
  self_loops = FALSE, q_max = Inf)

latent_multigraph_block_state(graph, nested = TRUE, blocks = NULL,
  latent_weights = NULL, degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"),
  self_loops = TRUE, q_max = Inf)

uncertain_description_length(state, sbm = TRUE, density = FALSE,
  aE = 1, breakdown = FALSE)

measured_error_posterior(state)

edge_probabilities(state, edges, log = TRUE, epsilon = 1e-8)

reconstruct_network(state, nsamples = 200L, burnin = 100L, thin = 1L,
  beta = 1, pedges = 0.5, ph = 0.1, hstep = 1,
  c = 1, d = 0.01, init_sbm = TRUE,
  collect_partitions = FALSE, threshold = 0.5)

latent_multigraph(graph, epsilon = 1e-8, max_niter = 0L)
}
\arguments{
  \item{graph}{A graph created by \code{\link{graph_create}}.}
  \item{n, x}{Total measurement trials \eqn{n_{ij}} and positive observations
  \eqn{x_{ij}} for each observed edge of \code{graph}, supplied as a numeric
  vector of length \code{graph_num_edges(graph)}, an edge property name, or
  \code{NULL} (defaulting to 1).}
  \item{n_default, x_default}{Default measurement counts for unobserved dyads
  not listed in \code{graph}.}
  \item{fn_params}{Named list \code{list(alpha, beta)} giving Beta prior
  parameters for the false-negative error rate \eqn{p}.}
  \item{fp_params}{Named list \code{list(mu, nu)} giving Beta prior parameters
  for the false-positive error rate \eqn{q}.}
  \item{lp, lq}{Optional fixed log error rates \eqn{\ln p} and \eqn{\ln q}
  (\eqn{\le 0}). Leave as \code{NA_real_} to integrate out unknown error rates
  under their Beta priors.}
  \item{q}{Extraneous edge probability estimates \eqn{q_{ij} \in [0, 1]} for
  observed edges of \code{graph}, given as a numeric vector, edge property
  name, or \code{NULL}.}
  \item{q_default}{Default probability in \eqn{[0, 1]} for unobserved dyads.}
  \item{nested}{Whether to place a nested hierarchical SBM prior (\code{TRUE})
  or a flat SBM prior (\code{FALSE}) on the latent multigraph.}
  \item{blocks}{Optional initial partition: \code{NULL}, a length-\eqn{N}
  integer vector, a hierarchy list, or a fitted block state.}
  \item{latent_edges, latent_weights}{Optional initial latent multigraph
  edgelist (2-column matrix) and integer multiplicities. Defaults to the edges
  of \code{graph} with unit weight.}
  \item{degree_corrected}{Whether the SBM prior is degree-corrected.}
  \item{degree_prior}{Degree sequence hyperprior (\code{"distributed"} or
  \code{"uniform"}).}
  \item{self_loops}{Whether self-loops are valid dyads in the latent graph.}
  \item{q_max}{Maximum degree cached in the partition-count table.}
  \item{state}{A \code{rgraphtool_uncertain_block_state} object.}
  \item{sbm}{Include the SBM prior description length (\code{TRUE}) or return
  only the measurement negative log-likelihood (\code{FALSE}).}
  \item{density}{Include an optional Poisson edge-count density prior with mean
  \code{aE}.}
  \item{aE}{Expected total number of edges when \code{density = TRUE}.}
  \item{breakdown}{If \code{TRUE}, return a named list separating total,
  measurement, and SBM description lengths plus Beta posterior parameters.}
  \item{edges}{Two-column matrix or data.frame of 1-based vertex pairs to query.}
  \item{log}{Return log-probabilities (\code{TRUE}) or probabilities
  (\code{FALSE}).}
  \item{epsilon}{Convergence tolerance for latent multiplicity summation or EM
  iteration.}
  \item{nsamples, burnin, thin}{Number of posterior samples collected, burn-in
  sweeps discarded, and sweeps between consecutive samples.}
  \item{beta}{Inverse MCMC temperature (default 1 for posterior sampling).}
  \item{pedges}{Probability of performing an edge-multiplicity sweep versus an
  SBM partition sweep.}
  \item{ph, hstep}{Hyperparameter move probability and proposal step size for
  \code{mixed_measured_block_state()}.}
  \item{c, d}{SBM proposal mixing and empty-group creation parameters.}
  \item{init_sbm}{If \code{TRUE} and \code{state} starts from a single-group
  partition, initialize the SBM hierarchy on the initial latent graph before
  MCMC.}
  \item{collect_partitions}{If \code{TRUE}, return sampled SBM partitions in
  \code{partitions}.}
  \item{threshold}{Marginal edge-probability threshold used to construct the
  returned \code{consensus_graph}.}
  \item{max_niter}{Maximum EM iterations for \code{latent_multigraph()}
  (\code{0L} means run until \code{epsilon} convergence).}
}
\description{
Native C++ implementations of Bayesian network reconstruction from noisy or
uncertain measurements (Peixoto 2018) and erased Poisson multigraph inference
(Peixoto 2020), matching graph-tool's \code{MeasuredBlockState},
\code{MixedMeasuredBlockState}, \code{UncertainBlockState},
\code{LatentMultigraphBlockState}, and \code{latent_multigraph()}.
}
\details{
All four reconstruction models combine an exact dynamic-edge flat or nested SBM
prior with a measurement likelihood over observed dyads:
\itemize{
  \item \code{measured_block_state()}: repeated binary observations
  \eqn{(n_{ij}, x_{ij})} with unknown Beta-prior or fixed homogeneous
  false-negative (\eqn{p}) and false-positive (\eqn{q}) rates.
  \item \code{mixed_measured_block_state()}: heterogeneous per-dyad
  Beta-Binomial error rates with joint hyperparameter MCMC over
  \eqn{(\alpha, \beta, \mu, \nu)}.
  \item \code{uncertain_block_state()}: extraneous per-edge probability
  estimates \eqn{q_{ij}} and default non-edge probability \eqn{q_{\mathrm{default}}}.
  \item \code{latent_multigraph_block_state()}: Bayesian posterior over latent
  edge multiplicities and self-loops conditioned on an erased simple graph.
}
\code{latent_multigraph()} implements the fast \eqn{O(V + E)}
Expectation-Maximization algorithm of Peixoto (2020) for maximum-likelihood
expected edge and self-loop multiplicities.
}
\value{
Constructors return a \code{rgraphtool_uncertain_block_state}.
\code{uncertain_description_length()} returns a numeric scalar (or breakdown
list). \code{measured_error_posterior()} returns Beta posterior parameters,
means, and standard deviations for \eqn{p} and \eqn{q}.
\code{edge_probabilities()} returns a numeric vector of conditional (log-)
probabilities. \code{reconstruct_network()} returns a
\code{rgraphtool_reconstruction} containing posterior marginal edge
probabilities (\code{edges}), a thresholded \code{consensus_graph}, the final
\code{state}, posterior error-rate traces, and entropy history.
\code{latent_multigraph()} returns a \code{rgraphtool_latent_multigraph} with
expected edge and self-loop weights (\code{edges}), vertex propensities
\code{theta_out}/\code{theta_in}, and convergence diagnostics.
}
\references{
Peixoto, T. P. (2018). Reconstructing networks with unknown and heterogeneous
errors. \emph{Physical Review X} 8, 041011. doi:10.1103/PhysRevX.8.041011

Peixoto, T. P. (2020). Latent Poisson models for networks with heterogeneous
density. \emph{Physical Review E} 102, 012309. doi:10.1103/PhysRevE.102.012309
}
\seealso{\code{\link{predict_edges}}, \code{\link{evaluate_edge_prediction}},
\code{\link{latent_closure_block_state}},
\code{demo("uncertain_reconstruction")}}
\examples{
g <- graph_create(8, directed = FALSE)
e <- rbind(t(combn(1:4, 2)), t(combn(5:8, 2)), c(4, 5))
graph_add_edges(g, e, properties = list(n = rep(3, nrow(e)), x = rep(3, nrow(e))))
st <- measured_block_state(g, n = "n", x = "x", n_default = 3L, x_default = 0L)
uncertain_description_length(st)
em <- latent_multigraph(g)
head(em$edges)
}
