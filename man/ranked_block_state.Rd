\name{ranked_block_state}
\alias{ranked_block_state}
\alias{ranked_description_length}
\alias{minimize_ranked_blockmodel_dl}
\alias{nested_ranked_block_state}
\alias{nested_ranked_level_entropy}
\alias{nested_ranked_description_length}
\alias{minimize_nested_ranked_blockmodel_dl}
\alias{ranked_vertex_order}
\alias{ranked_edge_directions}
\alias{ranked_edge_counts}
\alias{compare_ordered_models}
\alias{print.graphtool_ranked_block_state}
\alias{print.graphtool_nested_ranked_block_state}
\alias{print.graphtool_ordered_model_comparison}
\title{Ordered (ranked) flat and nested stochastic block models}
\usage{
ranked_block_state(graph, blocks = NULL, u = NULL, distance = FALSE,
  unique = FALSE, M = NULL, degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"))
ranked_description_length(state, label_dl = TRUE, degree_prior = NULL)
minimize_ranked_blockmodel_dl(graph, distance = FALSE, unique = FALSE,
  M = NULL, degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"), B_min = 1L, B_max = Inf,
  refine = TRUE, epsilon = 1e-3)
nested_ranked_block_state(graph, levels = NULL, u = NULL, distance = FALSE,
  unique = FALSE, M = NULL, degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"))
nested_ranked_level_entropy(state, degree_prior = NULL)
nested_ranked_description_length(state, degree_prior = NULL)
minimize_nested_ranked_blockmodel_dl(graph, distance = FALSE, unique = FALSE,
  M = NULL, degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"), B_min = 1L, B_max = Inf,
  refine = TRUE, epsilon = 1e-3)
ranked_vertex_order(state, normalize = FALSE)
ranked_edge_directions(state)
ranked_edge_counts(state)
compare_ordered_models(graph,
  models = c("unordered_flat", "unordered_nested",
    "ranked_flat", "ranked_nested"),
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  method = c("map", "posterior"),
  evidence_method = c("variational", "importance"),
  nruns = 1L, nsamples = 200L, thin = 5L, burnin = NULL, mle = TRUE,
  B_min = 1L, B_max = Inf, refine = TRUE)
}
\arguments{
  \item{graph}{A directed graph created by \code{\link{graph_create}} with
  \code{directed = TRUE}.}
  \item{blocks}{Optional integer vector of 1-based group memberships of length
  \code{graph_num_vertices(graph)}. Defaults to a single group.}
  \item{levels}{Optional list of at least two 1-based integer group membership
  vectors defining a nested hierarchy over the directed graph.}
  \item{u}{Optional finite integer vector of length at least \eqn{B} assigning
  each level-0 group an ordinal rank. Defaults to \code{0:(B - 1)}.}
  \item{distance}{If \code{FALSE} (the default), edge orientation probabilities
  depend only on whether \eqn{u_r < u_s}, \eqn{u_r = u_s}, or \eqn{u_r > u_s};
  if \code{TRUE}, edge orientation counts are stratified by rank distance
  \eqn{|u_r - u_s|}.}
  \item{unique}{If \code{FALSE} (the default), multiple groups may share the
  same rank drawn from \eqn{\{0, \dots, M - 1\}}; if \code{TRUE}, all groups
  have strictly distinct ranks.}
  \item{M}{Size of the discrete rank alphabet when \code{unique = FALSE}.
  Defaults to the number of vertices \eqn{N}.}
  \item{degree_corrected}{Whether to fit the degree-corrected (\code{TRUE}) or
  non-degree-corrected (\code{FALSE}) ranked SBM.}
  \item{degree_prior}{Degree sequence hyperprior (\code{"distributed"} or
  \code{"uniform"}).}
  \item{state}{A flat \code{graphtool_ranked_block_state} or nested
  \code{graphtool_nested_ranked_block_state} object.}
  \item{label_dl}{Whether to include the rank-label prior term
  \eqn{\ln \binom{M + B - 1}{B}} (when \code{unique = FALSE}).}
  \item{B_min, B_max}{Minimum and maximum number of level-0 occupied groups.}
  \item{refine}{Whether to run greedy refinement and rank-optimization sweeps
  after merge-split minimization.}
  \item{epsilon}{Stopping tolerance in nats for greedy refinement.}
  \item{normalize}{If \code{TRUE}, rescale vertex ranks to \eqn{[0, 1]}.}
  \item{models}{Character vector naming a subset of \code{"unordered_flat"},
  \code{"unordered_nested"}, \code{"ranked_flat"}, \code{"ranked_nested"},
  \code{"ranked_distance_flat"}, and \code{"ranked_distance_nested"}.}
  \item{method}{Compare models by \code{"map"} minimum description length or
  full Bayesian \code{"posterior"} marginal evidence.}
  \item{evidence_method}{Evidence estimator forwarded to
  \code{\link{model_evidence}} when \code{method = "posterior"}.}
  \item{nruns}{Number of independent minimization restarts per model.}
  \item{nsamples, thin, burnin, mle}{Posterior sampling and entropy estimation
  controls forwarded to \code{\link{model_evidence}}.}
}
\description{
Native C++ implementation of Peixoto's (2022) flat and nested ordered (ranked)
stochastic block models for directed networks, matching graph-tool's
\code{RankedBlockState} and
\code{NestedBlockState(base_state = RankedBlockState)}. Jointly infers
community assignments and an ordinal ranking of groups that explains global
directional flow (upstream, lateral, and downstream edges).
}
\details{
In the ranked SBM, the undirected multigraph of edge counts
\eqn{m_{rs} = e_{rs} + e_{sr}} is generated by an undirected SBM (or an
undirected nested hierarchy at levels \eqn{l \ge 1}), while edge orientations
between groups \eqn{r} and \eqn{s} are governed by their relative ranks
\eqn{u_r} and \eqn{u_s}. When \code{distance = FALSE}, minimization
canonicalizes distinct group ranks into contiguous integers \eqn{0, \dots, K-1}
oriented so that upstream edges (\eqn{u_r < u_s}) outnumber downstream edges
(\eqn{u_r > u_s}) without altering the exact description length.

Flat and nested ranked states integrate with \code{\link{block_membership}},
\code{\link{block_sizes}}, \code{\link{nested_levels}},
\code{\link{nested_project}}, \code{\link{nested_summary}},
\code{\link{mcmc_sweep}}, \code{\link{mcmc_equilibrate}},
\code{\link{mcmc_anneal}}, \code{\link{sample_posterior}},
\code{\link{model_evidence}}, \code{\link{plot_ranked_hierarchy}}, and
\code{\link{plot_hierarchy}}.
}
\value{
\code{ranked_block_state()} and \code{minimize_ranked_blockmodel_dl()} return a
\code{graphtool_ranked_block_state} list. \code{nested_ranked_block_state()} and
\code{minimize_nested_ranked_blockmodel_dl()} return a
\code{graphtool_nested_ranked_block_state} list.

\code{ranked_description_length()} and \code{nested_ranked_description_length()}
return total description lengths in nats; \code{nested_ranked_level_entropy()}
returns per-level contributions.

\code{ranked_vertex_order()} returns a numeric vector of per-vertex ranks.
\code{ranked_edge_directions()} returns an integer vector with \code{1L}
(upstream, \eqn{u_i < u_j}), \code{0L} (lateral, \eqn{u_i = u_j}), and
\code{-1L} (downstream, \eqn{u_i > u_j}) for every edge.
\code{ranked_edge_counts()} returns a list with counts \code{upstream},
\code{lateral}, \code{downstream}, and \code{internal}.

\code{compare_ordered_models()} returns a
\code{graphtool_ordered_model_comparison} list with ranked \code{summary},
\code{best_model}, \code{best_state}, \code{states}, optional \code{evidences},
and \code{log_bayes_factors}.
}
\references{
Peixoto, T. P. (2022). Ordered community structure in networks.
\emph{Physical Review E} 106, 024305. doi:10.1103/PhysRevE.106.024305
}
\seealso{\code{\link{plot_ranked_hierarchy}}, \code{\link{plot_hierarchy}},
\code{\link{mcmc_sweep}}, \code{\link{model_evidence}},
\code{demo("ordered_sbm")}}
\examples{
set.seed(1)
tier <- rep(1:3, each = 5L)
ed <- list()
for (i in 1:15) for (j in 1:15) if (i != j) {
  p <- if (tier[i] == tier[j]) 0.7 else if (tier[j] > tier[i]) 0.35 else 0.02
  if (runif(1) < p) ed[[length(ed) + 1L]] <- c(i, j)
}
g <- graph_add_edges(graph_create(15, directed = TRUE), do.call(rbind, ed))
st <- minimize_ranked_blockmodel_dl(g)
st
ranked_edge_counts(st)
}
