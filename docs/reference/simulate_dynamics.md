# Native C++ forward simulators for network dynamics

Fast native C++ simulators for Ising, compartmental epidemic
(SI/SIS/SIR/SIRS and SEI/SEIR), Gaussian, and generalized Potts dynamics
on networks, matching the generative likelihoods of
[`pseudo_ising_block_state`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md),
[`epidemics_block_state`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md),
[`pseudo_normal_block_state`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md),
[`linear_normal_block_state`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md),
and
[`pseudo_potts_block_state`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md).

## Usage

``` r
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
```

## Arguments

- graph:

  A `rgraphtool_graph` on which to simulate dynamics.

- weights:

  Edge coupling weight(s): scalar, edge property name, or numeric vector
  of length `graph_num_edges(graph)`.

- theta:

  Node external field(s) (scalar, length-\\N\\ vector, or \\N \times q\\
  matrix for Potts).

- n_samples, n_steps:

  Number of snapshots or time steps to record.

- burnin, thin:

  Burn-in sweeps discarded and sweeps between consecutive recorded
  snapshots.

- mode:

  Update schedule: `"equilibrium"` (asynchronous Gibbs sweeps),
  `"glauber"` (synchronous discrete-time updates), or `"linear"` (linear
  Gaussian difference equation).

- has_zero:

  For `simulate_ising()`, sample spins in `{-1, 0, +1}` instead of
  `{-1, +1}`.

- init_state, init_infected:

  Optional initial vertex states (or vector of initially infected
  1-based vertex IDs for epidemics).

- transmission:

  Edge transmission probability in `[0, 1)`.

- spontaneous:

  Per-step spontaneous infection probability in `[0, 1)`.

- recovery, resusceptible:

  Per-step recovery (`1 -> 2` or `1 -> 0` in SIS) and waning-immunity
  (`2 -> 0` in SIRS) probabilities.

- model:

  Compartmental epidemic model (`"SI"`, `"SIS"`, `"SIR"`, or `"SIRS"`).

- n_cascades:

  Number of independent epidemic cascades to simulate; if `> 1`, returns
  a list of \\N \times \mathrm{n\\steps}\\ matrices.

- exposed, expose_rate:

  Include an exposed compartment (`-1`) that transitions to infectious
  (`1`) with probability `expose_rate`.

- sigma, mu:

  Node standard deviation(s) and mean(s) for Gaussian dynamics.

- f:

  Square \\q \times q\\ interaction matrix for Potts dynamics.

## Value

An \\N \times M\\ integer or numeric matrix of simulated node states (or
a list of `n_cascades` matrices when `simulate_epidemics()` is called
with `n_cascades > 1`).

## See also

[`pseudo_ising_block_state`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md),
[`reconstruct_dynamics`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md),
`demo("dynamics_reconstruction")`

## Examples

``` r
set.seed(1)
g <- graph_create(6, directed = FALSE)
graph_add_edges(g, rbind(cbind(1:5, 2:6), c(6L, 1L)))
s <- simulate_ising(g, weights = 0.5, n_samples = 20L, burnin = 20L)
dim(s)
#> [1]  6 20
```
