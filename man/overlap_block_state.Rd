\name{overlap_block_state}
\alias{overlap_block_state}
\alias{overlap_description_length}
\alias{overlap_groups}
\alias{overlap_majority_groups}
\alias{minimize_overlap_blockmodel_dl}
\alias{nested_overlap_block_state}
\alias{nested_overlap_level_entropy}
\alias{nested_overlap_description_length}
\alias{minimize_nested_overlap_blockmodel_dl}
\alias{overlap_mcmc_sweep}
\alias{print.graphtool_overlap_block_state}
\alias{print.graphtool_nested_overlap_block_state}
\title{Flat and hierarchical overlapping stochastic block models}
\usage{
overlap_block_state(graph, half_edges = NULL, degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"), gen_overlap = TRUE)
overlap_description_length(state, degree_prior = NULL)
overlap_groups(state)
overlap_majority_groups(state)
minimize_overlap_blockmodel_dl(graph, degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"), gen_overlap = TRUE,
  B_min = 1L, B_max = Inf, refine = TRUE, epsilon = 1e-3)
nested_overlap_block_state(graph, half_edges, upper_levels,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  gen_overlap = TRUE)
nested_overlap_level_entropy(state, degree_prior = NULL)
nested_overlap_description_length(state, degree_prior = NULL)
minimize_nested_overlap_blockmodel_dl(graph, degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"), gen_overlap = TRUE,
  B_min = 1L, B_max = Inf, refine = TRUE, epsilon = 1e-3)
overlap_mcmc_sweep(state, beta = 1, niter = 1L, c = 0.5, d = 0.01,
  thin = 1L, collect = FALSE, B_min = 1L, B_max = Inf)
}
\arguments{
  \item{graph}{A graph created by \code{\link{graph_create}}.}
  \item{half_edges}{An \code{E x 2} integer matrix of positive group labels for
  the source and target half-edges of each row in \code{\link{graph_edges}(graph)},
  or \code{NULL} to place all half-edges in a single group.}
  \item{upper_levels}{Optional list of 1-based integer label vectors defining
  hierarchical partitions above the overlapping base level.}
  \item{degree_corrected}{Logical indicating whether to use the degree-corrected
  overlapping SBM.}
  \item{degree_prior}{Degree sequence hyperprior: \code{"distributed"} (default,
  Peixoto's hyperprior) or \code{"uniform"}.}
  \item{gen_overlap}{Logical indicating whether to use the generative overlapping
  mixture and degree hyperpriors (default \code{TRUE}, matching graph-tool 3.9).}
  \item{state}{An object returned by \code{overlap_block_state()},
  \code{minimize_overlap_blockmodel_dl()}, \code{nested_overlap_block_state()},
  or \code{minimize_nested_overlap_blockmodel_dl()}.}
  \item{B_min, B_max}{Minimum and maximum number of occupied level-0 groups.}
  \item{refine}{Whether to run greedy half-edge, vertex, edge-coupled, merge,
  split, and bridge-dissolution refinement sweeps.}
  \item{epsilon}{Relative entropy convergence tolerance for refinement.}
  \item{beta}{Inverse temperature for Metropolis-Hastings sweeps (\code{1} for
  exact Bayesian posterior sampling; \code{Inf} for greedy moves).}
  \item{niter}{Number of sweeps to perform.}
  \item{c, d}{Proposal mixing parameter and empty-group proposal probability.}
  \item{thin}{Keep every \code{thin}-th sweep when \code{collect = TRUE}.}
  \item{collect}{Whether to collect vertex majority-group partitions and
  description lengths during MCMC sweeps.}
}
\description{
Native C++ implementation of Peixoto's (2015) microcanonical overlapping
stochastic block model (\code{OverlapBlockState}) and its hierarchical extension
(\code{NestedBlockState} with an overlapping base level). Each edge endpoint
(half-edge) carries its own group label, allowing vertices to participate in
multiple communities simultaneously.
}
\details{
Exact flat and nested overlapping description lengths match graph-tool 3.9 to
machine precision for directed and undirected graphs, with or without parallel
edges, self-loops, isolated vertices, degree correction, and distributed or
uniform degree priors.

Minimization combines non-overlapping seeding, vertex- and edge-level
agglomeration, sequential bridge-group dissolution, and greedy half-edge,
vertex, and edge-coupled refinement sweeps. Per-vertex mixed memberships
returned by \code{overlap_groups(state)$matrix} can be plotted directly as
pie-chart nodes via \code{\link{plot_marginals}}.
}
\value{
\code{overlap_block_state()} and \code{minimize_overlap_blockmodel_dl()} return
a \code{graphtool_overlap_block_state} object containing \code{half_edges},
\code{overlap}, and model settings. \code{nested_overlap_block_state()} and
\code{minimize_nested_overlap_blockmodel_dl()} return a
\code{graphtool_nested_overlap_block_state} object adding \code{upper_levels}.
\code{overlap_groups()} returns a list with per-vertex \code{groups},
\code{out_degrees}, \code{in_degrees}, \code{total_degrees}, half-edge count
\code{matrix}, \code{majority} group vector, and occupied group count \code{B}.
\code{overlap_mcmc_sweep()} returns a list with updated \code{state},
\code{entropy}, \code{dS}, \code{nattempts}, \code{nmoves}, and (when
\code{collect = TRUE}) \code{majority_samples} and \code{entropies}.
}
\references{
Peixoto, T. P. (2015). Model selection and hypothesis testing for large-scale
network models with overlapping groups. \emph{Physical Review X} 5, 011033.
doi:10.1103/PhysRevX.5.011033
}
\seealso{\code{\link{compare_model_classes}}, \code{\link{model_evidence}},
\code{\link{plot_marginals}}, \code{demo("model_class_selection")}}
\examples{
g <- graph_create(n = 6L, directed = FALSE)
graph_add_edges(g, matrix(c(1, 2, 1, 3, 2, 3, 3, 4, 3, 5, 4, 5),
                          ncol = 2L, byrow = TRUE))
he <- matrix(c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2), ncol = 2L, byrow = TRUE)
st <- overlap_block_state(g, half_edges = he)
overlap_description_length(st)
overlap_groups(st)$matrix
}
