\name{dynamics_block_state}
\alias{pseudo_ising_block_state}
\alias{ising_glauber_block_state}
\alias{pseudo_cising_block_state}
\alias{cising_glauber_block_state}
\alias{epidemics_block_state}
\alias{pseudo_normal_block_state}
\alias{normal_glauber_block_state}
\alias{linear_normal_block_state}
\alias{pseudo_potts_block_state}
\alias{potts_glauber_block_state}
\alias{dynamics_description_length}
\alias{reconstruct_dynamics}
\alias{reconstruct_l1_path}
\alias{dynamics_precision_matrix}
\alias{dynamics_transmission_probs}
\alias{print.rgraphtool_dynamics_block_state}
\alias{print.rgraphtool_dynamics_reconstruction}
\alias{print.rgraphtool_l1_path}
\title{Network reconstruction from dynamics and behavior}
\usage{
pseudo_ising_block_state(s, graph = NULL, has_zero = FALSE,
  x = 1, x_range = c(-Inf, Inf), theta = 0, theta_range = c(-Inf, Inf),
  xdelta = 1e-8, tdelta = 1e-8, disable_xdist = FALSE, disable_tdist = FALSE,
  xl1 = 1, tl1 = 1, sbm = NULL, nested = TRUE, blocks = NULL,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  self_loops = FALSE, q_max = Inf)

ising_glauber_block_state(s, graph = NULL, directed = TRUE, has_zero = FALSE,
  x = 1, x_range = c(-Inf, Inf), theta = 0, theta_range = c(-Inf, Inf),
  xdelta = 1e-8, tdelta = 1e-8, disable_xdist = FALSE, disable_tdist = FALSE,
  xl1 = 1, tl1 = 1, sbm = NULL, nested = TRUE, blocks = NULL,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  self_loops = FALSE, q_max = Inf)

pseudo_cising_block_state(s, graph = NULL,
  x = 1, x_range = c(-Inf, Inf), theta = 0, theta_range = c(-Inf, Inf),
  xdelta = 1e-8, tdelta = 1e-8, disable_xdist = FALSE, disable_tdist = FALSE,
  xl1 = 1, tl1 = 1, sbm = NULL, nested = TRUE, blocks = NULL,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  self_loops = FALSE, q_max = Inf)

cising_glauber_block_state(s, graph = NULL, directed = TRUE,
  x = 1, x_range = c(-Inf, Inf), theta = 0, theta_range = c(-Inf, Inf),
  xdelta = 1e-8, tdelta = 1e-8, disable_xdist = FALSE, disable_tdist = FALSE,
  xl1 = 1, tl1 = 1, sbm = NULL, nested = TRUE, blocks = NULL,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  self_loops = FALSE, q_max = Inf)

epidemics_block_state(s, graph = NULL, directed = TRUE, exposed = FALSE,
  x = -0.5, x_range = c(-Inf, -1e-8), theta = -3, theta_range = c(-Inf, -1e-8),
  xdelta = 1e-8, tdelta = 1e-8, disable_xdist = FALSE, disable_tdist = FALSE,
  xl1 = 1, tl1 = 1, sbm = NULL, nested = TRUE, blocks = NULL,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  self_loops = FALSE, q_max = Inf)

pseudo_normal_block_state(s = NULL, graph = NULL, fix_mean = FALSE,
  positive = TRUE, pslack = 1e-6,
  x = 0.1, x_range = c(-Inf, Inf), theta = c(0, 0),
  theta_range = list(c(-200, 200), c(-Inf, Inf)),
  xdelta = 1e-8, tdelta = c(1e-8, 1e-8),
  disable_xdist = FALSE, disable_tdist = FALSE,
  xl1 = 1, tl1 = c(1, 1), sbm = NULL, nested = TRUE, blocks = NULL,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  self_loops = FALSE, q_max = Inf, cov = NULL, n_eff = NULL, mean = NULL)

normal_glauber_block_state(s, graph = NULL, directed = TRUE, fix_mean = FALSE,
  x = 0.1, x_range = c(-Inf, Inf), theta = c(0, 0),
  theta_range = list(c(-200, 200), c(-Inf, Inf)),
  xdelta = 1e-8, tdelta = c(1e-8, 1e-8),
  disable_xdist = FALSE, disable_tdist = FALSE,
  xl1 = 1, tl1 = c(1, 1), sbm = NULL, nested = TRUE, blocks = NULL,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  self_loops = FALSE, q_max = Inf)

linear_normal_block_state(s, graph = NULL, directed = TRUE, self_loops = TRUE,
  x = 0.1, x_range = c(-Inf, Inf), theta = 0, theta_range = c(-200, 200),
  xdelta = 1e-8, tdelta = 1e-8, disable_xdist = FALSE, disable_tdist = FALSE,
  xl1 = 1, tl1 = 1, sbm = NULL, nested = TRUE, blocks = NULL,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  q_max = Inf)

pseudo_potts_block_state(s, f = NULL, q = NULL, graph = NULL,
  x = 1, x_range = c(-Inf, Inf), theta = NULL, theta_range = c(-Inf, Inf),
  xdelta = 1e-8, tdelta = 1e-8, disable_xdist = FALSE, disable_tdist = FALSE,
  xl1 = 1, tl1 = 1, sbm = NULL, nested = TRUE, blocks = NULL,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  self_loops = FALSE, q_max = Inf)

potts_glauber_block_state(s, f = NULL, q = NULL, graph = NULL, directed = TRUE,
  x = 1, x_range = c(-Inf, Inf), theta = NULL, theta_range = c(-Inf, Inf),
  xdelta = 1e-8, tdelta = 1e-8, disable_xdist = FALSE, disable_tdist = FALSE,
  xl1 = 1, tl1 = 1, sbm = NULL, nested = TRUE, blocks = NULL,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  self_loops = FALSE, q_max = Inf)

dynamics_description_length(state, sbm = NULL, latent_edges = TRUE,
  density = FALSE, aE = 1, xdist = NULL, tdist = NULL,
  xdist_uniform = FALSE, tdist_uniform = FALSE, xl1 = NULL, tl1 = NULL,
  normal = FALSE, mu = 0, sigma = 1, breakdown = FALSE)

reconstruct_dynamics(state, niter = 25L, nsamples = 0L, burnin = 50L,
  thin = 2L, beta = Inf, k_cand = 20L, w_step = 0.25, t_step = 0.25,
  tol = 1e-4, threshold = 0.5, sbm = NULL, xdist = NULL, tdist = NULL,
  xl1 = NULL, tl1 = NULL, normal = FALSE, mu = 0, sigma = 1,
  refresh_every = 3L, merge_split = TRUE)

reconstruct_l1_path(state, xl1_values = c(0.25, 0.5, 1, 2, 5, 10),
  tl1 = 0.1, sbm = FALSE, niter = 20L, truth = NULL, k_cand = 20L,
  tol = 1e-4)

dynamics_precision_matrix(state)

dynamics_transmission_probs(state)
}
\arguments{
  \item{s}{An \eqn{N \times M} numeric/integer matrix of node observations (or
  a list of \eqn{N \times M_c} matrices representing independent time-series
  cascades or restarts). See Details for per-model value and preprocessing
  requirements.}
  \item{graph}{Optional initial \code{rgraphtool_graph} with \eqn{N} vertices.
  If \code{NULL}, reconstruction starts from an empty graph.}
  \item{directed}{Whether the reconstructed graph is directed.}
  \item{has_zero}{For Ising models, allow three-state spins \code{-1, 0, 1}
  instead of binary \code{-1, 1}.}
  \item{exposed}{For epidemic models, use SEI/SEIR transitions where newly
  infected nodes transition from \code{0} (susceptible) to \code{-1} (exposed)
  before becoming \code{1} (infectious).}
  \item{fix_mean}{For Gaussian models, center each vertex's observations to
  zero mean across samples.}
  \item{cov, n_eff, mean}{For \code{pseudo_normal_block_state()}, sufficient
  statistics used instead of \code{s}: an \code{N x N} covariance matrix
  \code{cov}, the effective number of observations \code{n_eff}, and an
  optional mean vector (default zero, e.g. for regression residuals). Give
  exactly one of \code{s} or \code{cov}.}
  \item{positive, pslack}{For \code{pseudo_normal_block_state()}, enforce strict
  diagonal dominance (positive definiteness) of the reconstructed precision
  matrix with margin \code{pslack}.}
  \item{f, q}{For Potts models, the \eqn{q \times q} spin interaction matrix
  (defaulting to the \eqn{q \times q} identity) and number of discrete spin
  states \eqn{q \ge 2}.}
  \item{x, x_range}{Initial edge weight(s) (scalar, edge property name, or
  vector) and allowed weight interval \code{c(xmin, xmax)}.}
  \item{theta, theta_range}{Initial node parameters (scalar, length-\code{q_t}
  vector, or \eqn{N \times q_t} matrix) and allowed bounds.}
  \item{xdelta, tdelta}{Quantization precision parameters \eqn{\Delta} for edge
  weights and node parameters.}
  \item{disable_xdist, disable_tdist}{If \code{TRUE}, disable quantized MDL
  weight/node distributions and use continuous L1 (Laplace) regularization
  instead.}
  \item{xl1, tl1}{L1 / Laplace scale parameter(s) \eqn{\lambda} for edge
  weights and node parameters.}
  \item{sbm}{Include the flat or nested SBM structural prior (\code{TRUE}) or
  reconstruct without community structure (\code{FALSE}).}
  \item{nested, blocks, degree_corrected, degree_prior, self_loops, q_max}{SBM
  prior configuration options.}
  \item{state}{A \code{rgraphtool_dynamics_block_state} (or
  \code{rgraphtool_dynamics_reconstruction} for precision/transmission helpers).}
  \item{latent_edges, density, aE, xdist, tdist, xdist_uniform, tdist_uniform, normal, mu, sigma}{Prior
  and description-length component overrides matching graph-tool's
  \code{DynamicsBlockStateBase.entropy()}.}
  \item{breakdown}{If \code{TRUE}, return a named list separating \code{total},
  \code{nll}, \code{geom_prior}, \code{x_prior}, \code{t_prior}, and
  \code{sbm}.}
  \item{niter}{Number of greedy MAP coordinate-descent / candidate-edge sweeps.}
  \item{nsamples, burnin, thin, beta}{Optional posterior MCMC sampling controls
  (set \code{nsamples > 0} and finite \code{beta = 1} to collect posterior
  marginal edge probabilities and weight/parameter uncertainties).}
  \item{k_cand}{Number of correlation-screened candidate partners per vertex
  for subquadratic scaling on larger networks.}
  \item{w_step, t_step}{Random-walk proposal half-widths for finite-temperature
  edge weight and node parameter MCMC updates.}
  \item{tol}{Stopping tolerance in nats for greedy MAP optimization.}
  \item{refresh_every}{Active-set candidate refresh period (in greedy MAP iterations) using residual correlations, 2-hop neighbors, and SBM block proposals; set \code{0L} to keep initial candidates fixed.}
  \item{merge_split}{If \code{TRUE} and \code{sbm = TRUE}, perform multilevel merge-split sweeps inside \code{DynBlockPrior} during greedy MAP optimization.}
  \item{threshold}{Marginal edge-probability threshold for the returned
  \code{graph}.}
  \item{xl1_values}{Grid of L1 regularization strengths \eqn{\lambda} swept
  with warm starts in \code{reconstruct_l1_path()}.}
  \item{truth}{Optional ground-truth \code{rgraphtool_graph} or 2-column
  edgelist matrix used to compute precision, recall, F1, and ROC-AUC along the
  L1 path.}
}
\description{
Native C++ implementations of indirect network reconstruction from dynamical
time series and equilibrium behavioral snapshots (Peixoto 2019, 2024, 2025),
supporting all 10 graph-tool dynamical/behavioral models, L1 (Laplace),
Gaussian, and quantized MDL weight/node priors, joint flat/nested SBM community
detection, and declarative L1 regularization path sweeps.
}
\details{
Input and preprocessing conventions by model family:
\itemize{
  \item \strong{Ising} (\code{pseudo_ising_block_state},
  \code{ising_glauber_block_state}): integer spins in \code{\{-1, +1\}} (or
  \code{\{-1, 0, +1\}} with \code{has_zero = TRUE}); continuous Ising models
  accept real values in \code{[-1, 1]}.
  \item \strong{Epidemics} (\code{epidemics_block_state}): integer states with
  \code{0} = susceptible, \code{1} = infectious (and \code{-1} = exposed when
  \code{exposed = TRUE}; recovered states such as \code{2} are handled
  automatically). Edge weights \eqn{x_{ji} \le 0} encode transmission
  probabilities \eqn{t_{ji} = 1 - e^{x_{ji}}} and node parameters
  \eqn{\theta_i \le 0} encode spontaneous infection probabilities
  \eqn{r_i = e^{\theta_i}}.
  \item \strong{Continuous Gaussian} (\code{pseudo_normal_block_state},
  \code{normal_glauber_block_state}, \code{linear_normal_block_state}): real
  observations. For empirical financial or sensor time series, convert levels
  to stationary log-returns or first differences and standardize each variable
  to zero mean and unit variance before fitting.
  \item \strong{Potts} (\code{pseudo_potts_block_state},
  \code{potts_glauber_block_state}): integer categories in \code{0..(q - 1)}.
}
}
\value{
Constructors return a \code{rgraphtool_dynamics_block_state}.
\code{dynamics_description_length()} returns a scalar DL (or component
breakdown list). \code{reconstruct_dynamics()} returns a
\code{rgraphtool_dynamics_reconstruction} containing the inferred \code{edges}
data.frame, reconstructed \code{graph}, estimated \code{theta} and
\code{theta_sd}, updated \code{state}, and entropy trace.
\code{reconstruct_l1_path()} returns a \code{rgraphtool_l1_path} with ranked
\code{summary}, per-lambda \code{fits}, and BIC-selected \code{best_fit}.
\code{dynamics_precision_matrix()} returns the precision matrix, partial
correlation matrix, \code{sigma}, and \code{mu}.
\code{dynamics_transmission_probs()} returns edge transmission probabilities
and per-node spontaneous infection probabilities.
}
\references{
Peixoto, T. P. (2019). Network reconstruction and community detection from
dynamics. \emph{Physical Review Letters} 123, 128301.
doi:10.1103/PhysRevLett.123.128301

Peixoto, T. P. (2025). Network reconstruction via the minimum description
length principle. \emph{Physical Review X} 15, 011065.
doi:10.1103/PhysRevX.15.011065
}
\seealso{\code{\link{simulate_ising}}, \code{\link{simulate_epidemics}},
\code{\link{simulate_normal_dynamics}}, \code{\link{simulate_potts}},
\code{\link{evaluate_edge_prediction}}, \code{demo("dynamics_reconstruction")},
\code{demo("empirical_dynamics")}}
\examples{
set.seed(1)
g <- graph_create(8, directed = FALSE)
graph_add_edges(g, rbind(cbind(1:7, 2:8), c(8L, 1L)))
s <- simulate_ising(g, weights = 0.6, n_samples = 200L, burnin = 100L, thin = 2L)
st <- pseudo_ising_block_state(s, nested = FALSE, xdelta = 1e-4, tdelta = 1e-4)
rec <- reconstruct_dynamics(st, niter = 10L, xdist = FALSE, tdist = FALSE, sbm = TRUE)
rec
}
