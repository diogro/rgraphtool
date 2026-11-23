\name{simulate_dynamics}
\alias{simulate_ising}
\alias{simulate_epidemics}
\alias{simulate_normal_dynamics}
\alias{simulate_potts}
\title{Native C++ forward simulators for network dynamics}
\usage{
simulate_ising(graph, weights = 0.5, theta = 0, n_samples = 200L,
  burnin = 100L, thin = 1L, mode = c("equilibrium", "glauber"),
  has_zero = FALSE, init_state = NULL)

simulate_epidemics(graph, transmission = 0.3, spontaneous = 0.005,
  recovery = 0.2, resusceptible = 0.0,
  model = c("SI", "SIS", "SIR", "SIRS"),
  n_steps = 50L, n_cascades = 1L, exposed = FALSE, expose_rate = 0.5,
  init_infected = NULL)

simulate_normal_dynamics(graph, weights = 0.1, sigma = 1, mu = 0,
  n_samples = 200L, burnin = 100L, thin = 1L,
  mode = c("equilibrium", "glauber", "linear"))

simulate_potts(graph, weights = 0.5, f = diag(3), theta = NULL,
  n_samples = 200L, burnin = 100L, thin = 1L,
  mode = c("equilibrium", "glauber"), init_state = NULL)
}
\arguments{
  \item{graph}{A \code{rgraphtool_graph} on which to simulate dynamics.}
  \item{weights}{Edge coupling weight(s): scalar, edge property name, or
  numeric vector of length \code{graph_num_edges(graph)}.}
  \item{theta}{Node external field(s) (scalar, length-\eqn{N} vector, or
  \eqn{N \times q} matrix for Potts).}
  \item{n_samples, n_steps}{Number of snapshots or time steps to record.}
  \item{burnin, thin}{Burn-in sweeps discarded and sweeps between consecutive
  recorded snapshots.}
  \item{mode}{Update schedule: \code{"equilibrium"} (asynchronous Gibbs
  sweeps), \code{"glauber"} (synchronous discrete-time updates), or
  \code{"linear"} (linear Gaussian difference equation).}
  \item{has_zero}{For \code{simulate_ising()}, sample spins in
  \code{\{-1, 0, +1\}} instead of \code{\{-1, +1\}}.}
  \item{init_state, init_infected}{Optional initial vertex states (or vector of
  initially infected 1-based vertex IDs for epidemics).}
  \item{transmission}{Edge transmission probability in \code{[0, 1)}.}
  \item{spontaneous}{Per-step spontaneous infection probability in
  \code{[0, 1)}.}
  \item{recovery, resusceptible}{Per-step recovery (\code{1 -> 2} or
  \code{1 -> 0} in SIS) and waning-immunity (\code{2 -> 0} in SIRS)
  probabilities.}
  \item{model}{Compartmental epidemic model (\code{"SI"}, \code{"SIS"},
  \code{"SIR"}, or \code{"SIRS"}).}
  \item{n_cascades}{Number of independent epidemic cascades to simulate; if
  \code{> 1}, returns a list of \eqn{N \times \mathrm{n\_steps}} matrices.}
  \item{exposed, expose_rate}{Include an exposed compartment (\code{-1}) that
  transitions to infectious (\code{1}) with probability \code{expose_rate}.}
  \item{sigma, mu}{Node standard deviation(s) and mean(s) for Gaussian
  dynamics.}
  \item{f}{Square \eqn{q \times q} interaction matrix for Potts dynamics.}
}
\description{
Fast native C++ simulators for Ising, compartmental epidemic (SI/SIS/SIR/SIRS
and SEI/SEIR), Gaussian, and generalized Potts dynamics on networks, matching
the generative likelihoods of \code{\link{pseudo_ising_block_state}},
\code{\link{epidemics_block_state}}, \code{\link{pseudo_normal_block_state}},
\code{\link{linear_normal_block_state}}, and
\code{\link{pseudo_potts_block_state}}.
}
\value{
An \eqn{N \times M} integer or numeric matrix of simulated node states (or a
list of \code{n_cascades} matrices when \code{simulate_epidemics()} is called
with \code{n_cascades > 1}).
}
\seealso{\code{\link{pseudo_ising_block_state}},
\code{\link{reconstruct_dynamics}}, \code{demo("dynamics_reconstruction")}}
\examples{
set.seed(1)
g <- graph_create(6, directed = FALSE)
graph_add_edges(g, rbind(cbind(1:5, 2:6), c(6L, 1L)))
s <- simulate_ising(g, weights = 0.5, n_samples = 20L, burnin = 20L)
dim(s)
}
