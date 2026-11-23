\name{latent_closure_block_state}
\alias{latent_closure_block_state}
\alias{closure_description_length}
\alias{infer_triadic_closure}
\alias{print.rgraphtool_latent_closure_block_state}
\alias{print.rgraphtool_closure_inference}
\title{Latent triadic closure stochastic block model}
\usage{
latent_closure_block_state(graph, L = 1L, nested = TRUE,
  blocks = NULL, gen = NULL, w0 = NULL, loops = NULL, ex = NULL,
  degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"), q_max = Inf)

closure_description_length(state, breakdown = FALSE)

infer_triadic_closure(state, nsamples = 200L, burnin = 100L, thin = 1L,
  beta = 1, pedges = 0.5, c = 1, d = 0.01,
  sample_loops = TRUE, init_sbm = TRUE,
  collect_partitions = FALSE)
}
\arguments{
  \item{graph}{An undirected simple graph created by \code{\link{graph_create}}
  (no self-loops or parallel edges).}
  \item{L}{Maximum number of triadic closure generations (\eqn{L \ge 1}).}
  \item{nested}{Whether to use a nested (\code{TRUE}) or flat (\code{FALSE})
  SBM substrate prior.}
  \item{blocks}{Optional initial substrate partition or hierarchy.}
  \item{gen}{Optional integer vector of length \code{graph_num_edges(graph)} in
  \code{0..L} giving each edge's initial generation (\code{0} = substrate,
  \code{1..L} = triadic closure generation). Defaults to all \code{0L}.}
  \item{w0}{Optional initial substrate edge multiplicities for generation-0
  edges.}
  \item{loops}{Optional length-\eqn{N} integer vector of initial substrate
  self-loop multiplicities.}
  \item{ex}{Optional list of length \code{graph_num_edges(graph)} giving the
  1-based closing ego vertex IDs for each closure edge (empty integer vectors
  for substrate edges).}
  \item{degree_corrected}{Whether the substrate SBM is degree-corrected.}
  \item{degree_prior}{Degree sequence hyperprior (\code{"distributed"} or
  \code{"uniform"}).}
  \item{q_max}{Maximum degree cached in the partition-count table.}
  \item{state}{A \code{rgraphtool_latent_closure_block_state} object.}
  \item{breakdown}{If \code{TRUE}, return a named list with \code{total},
  \code{sbm}, \code{closure}, per-generation \code{layer_entropies}, and
  per-vertex open/closed triad counts \code{M} and \code{E}.}
  \item{nsamples, burnin, thin}{Posterior MCMC sample count, burn-in sweeps,
  and thinning interval.}
  \item{beta}{Inverse MCMC temperature.}
  \item{pedges}{Probability of performing an edge-layer/ego sweep versus a
  substrate SBM partition sweep.}
  \item{c, d}{Substrate SBM proposal parameters.}
  \item{sample_loops}{Whether to sample latent substrate self-loop
  multiplicities during edge sweeps.}
  \item{init_sbm}{If \code{TRUE} and \code{state} starts from a single-group
  partition, initialize the substrate SBM partition before MCMC.}
  \item{collect_partitions}{If \code{TRUE}, return sampled substrate SBM
  partitions in \code{partitions}.}
}
\description{
Native C++ implementation of Peixoto's (2022) joint stochastic block model with
\eqn{L \ge 1} generations of zero-inflated triadic closures, matching
graph-tool's \code{LatentClosureBlockState}. Disentangles observed edges into
homophily/community substrate edges versus triadic-closure edges and estimates
per-vertex open/closed triad counts and transitivity.
}
\value{
\code{latent_closure_block_state()} returns a
\code{rgraphtool_latent_closure_block_state}.
\code{closure_description_length()} returns the exact joint description length
(or a detailed breakdown list). \code{infer_triadic_closure()} returns a
\code{rgraphtool_closure_inference} containing per-edge posterior layer
probabilities (\code{edges}: \code{substrate_prob}, \code{closure_prob}, and
\code{gen_0}..\code{gen_L}), per-generation vertex statistics
(\code{vertex_stats}: average closed triads \code{t}, open triads \code{m}, and
transitivity fraction \code{c}), the final \code{state}, and entropy trace.
}
\references{
Peixoto, T. P. (2022). Disentangling homophily, community structure, and
triadic closure in networks. \emph{Physical Review X} 12, 011004.
doi:10.1103/PhysRevX.12.011004
}
\seealso{\code{\link{measured_block_state}}, \code{\link{predict_edges}},
\code{demo("triadic_closure")}}
\examples{
g <- graph_create(6, directed = FALSE)
graph_add_edges(g, rbind(cbind(1:5, 2:6), c(6L, 1L), c(1L, 3L), c(3L, 5L)))
st <- latent_closure_block_state(g, L = 1L, nested = FALSE)
set.seed(1)
res <- infer_triadic_closure(st, nsamples = 40L, burnin = 30L, init_sbm = FALSE)
res$edges
}
