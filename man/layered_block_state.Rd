\name{layered_block_state}
\alias{layered_block_state}
\alias{layered_description_length}
\alias{minimize_layered_blockmodel_dl}
\alias{nested_layered_block_state}
\alias{nested_layered_level_entropy}
\alias{nested_layered_description_length}
\alias{minimize_nested_layered_blockmodel_dl}
\alias{print.rgraphtool_layered_block_state}
\alias{print.rgraphtool_nested_layered_block_state}
\title{Flat and nested layered stochastic block models}
\usage{
layered_block_state(graph, layers, independent = TRUE, blocks = NULL,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  weights = NULL, types = NULL, params = "microcanonical", present = NULL)
layered_description_length(state, lmemb = is.null(state$present),
  degree_prior = NULL)
minimize_layered_blockmodel_dl(graph, layers, independent = TRUE,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  B_min = 1L, B_max = Inf, refine = TRUE, epsilon = 1e-3,
  weights = NULL, types = NULL, params = "microcanonical",
  nthreads = NULL, present = NULL)
nested_layered_block_state(graph, layers, independent = TRUE, levels = NULL,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  weights = NULL, types = NULL, params = "microcanonical", present = NULL)
nested_layered_level_entropy(state, degree_prior = NULL)
nested_layered_description_length(state, degree_prior = NULL)
minimize_nested_layered_blockmodel_dl(graph, layers, independent = TRUE,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  B_min = 1L, B_max = Inf, refine = TRUE, epsilon = 1e-3,
  weights = NULL, types = NULL, params = "microcanonical",
  nthreads = NULL, present = NULL)
}
\arguments{
  \item{graph}{A graph created by \code{\link{graph_create}}.}
  \item{layers}{Discrete edge layer labels: either a character scalar naming an
  edge property stored on \code{graph} or a numeric, integer, character, factor,
  or logical vector of length \code{graph_num_edges(graph)}.}
  \item{independent}{If \code{TRUE} (default), fit the conditionally independent
  layers model (each layer has its own active vertices, degree sequence, and
  block edge counts conditioned on the shared partition). If \code{FALSE}, fit
  the edge-covariate model (a single collapsed topology and degree sequence
  with categorical layer labels per block pair).}
  \item{blocks}{Optional 1-based integer group membership vector of length
  \code{graph_num_vertices(graph)}.}
  \item{levels}{Optional non-empty list of 1-based integer label vectors for
  each hierarchy level.}
  \item{degree_corrected}{Whether to use degree correction (\code{TRUE}) or the
  non-degree-corrected model (\code{FALSE}).}
  \item{degree_prior}{Degree sequence hyperprior (\code{"distributed"} or
  \code{"uniform"}).}
  \item{state}{A \code{rgraphtool_layered_block_state} or
  \code{rgraphtool_nested_layered_block_state} object.}
  \item{lmemb}{Whether to include the layer-membership description length term
  (applicable when \code{independent = TRUE}). Defaults to \code{TRUE} when
  presence is inferred from the edges and \code{FALSE} when \code{present}
  was given, since given presence is not encoded.}
  \item{B_min, B_max}{Minimum and maximum number of occupied level-0 groups.}
  \item{refine}{Whether to run local refinement sweeps during minimization.}
  \item{epsilon}{Convergence tolerance in nats.}
  \item{weights}{Optional edge covariate values (same formats accepted by
  \code{\link{weighted_block_state}}) for combined layered-weighted SBMs.}
  \item{types}{Optional character vector of covariate distribution types
  (required when \code{weights} is supplied).}
  \item{params}{Optional covariate hyperprior specification (default
  \code{"microcanonical"}).}
  \item{nthreads}{Optional thread count override (defaults to
  \code{rgraphtool_num_threads()}).}
  \item{present}{Optional known per-layer vertex sets for the independent
  layers model: a logical matrix with one row per vertex and one column per
  layer (columns may be named by layer label). \code{NULL} (default) infers
  each layer's vertex set from its edges and encodes it (graph-tool's
  behaviour). When given, each layer is an SBM over its present vertices
  (which may have degree zero), presence is conditioned on rather than
  encoded, and the layer-membership term is dropped. Use it when vertex sets
  differ by design, e.g. genes measured in only some conditions. Every edge
  endpoint must be present in the edge's layer. Ignored when
  \code{independent = FALSE}.}
}
\description{
Construct, evaluate, and infer flat and hierarchical (nested) layered
stochastic block models (Peixoto 2015), supporting both the conditionally
independent layers model (\code{independent = TRUE}) and the discrete
edge-covariate model (\code{independent = FALSE}).
}
\details{
All description lengths match graph-tool 3.9's \code{LayeredBlockState} and
\code{NestedBlockState(base_state=LayeredBlockState)} to machine precision for
directed and undirected graphs (including parallel edges within and across
layers, self-loops, isolated vertices, and blocks absent from individual
layers).

With \code{present} given, independent-layers description lengths are exactly
normalised over graphs on the given vertex sets. Blocks with no present vertex
in a layer still reserve block-pair edge counts, so the code is slightly
conservative (sums below one) in that case.

With edge weights, the weight hyperprior of the independent layers model is
pooled over layers: every (layer, block pair) weight distribution is drawn from
one shared hyperprior, as the block pairs are in the edge-covariate model. Both
models then leave the same summary statistics unencoded, so their difference
does not depend on the units of the weights. This departs from graph-tool, which
gives every layer its own hyperprior and so shifts the independent-layers
description length by \eqn{(C - 1)\log s} (exponential weights, \eqn{C}
layers) when weights are multiplied by \eqn{s}. One offset remains for
\code{"real-normal"} weights when the top level holds a single block pair
(always in a nested hierarchy; a flat fit only at \eqn{B = 1}): the
edge-covariate model then has a single top-level weight sum, which has no spread
to leave unencoded, so the independent-minus-edge-covariate difference changes by
\eqn{-2\log s}, whatever the number of layers. Weights in natural, unit-free
scales (such as correlations) avoid the issue.

Layered states integrate directly with \code{\link{block_membership}},
\code{\link{block_sizes}}, \code{\link{nested_levels}},
\code{\link{nested_project}}, \code{\link{nested_summary}},
\code{\link{plot_hierarchy}}, \code{\link{mcmc_sweep}},
\code{\link{mcmc_equilibrate}}, \code{\link{mcmc_anneal}},
\code{\link{sample_posterior}}, and \code{\link{model_evidence}}.
}
\value{
Flat constructors and minimizers return a \code{rgraphtool_layered_block_state};
nested constructors and minimizers return a
\code{rgraphtool_nested_layered_block_state}. Description length functions return
a numeric scalar (or per-level numeric vector for
\code{nested_layered_level_entropy()}).
}
\references{
Peixoto, T. P. (2015). Inferring the mesoscale structure of layered,
edge-valued, and time-varying networks. \emph{Physical Review E} 92, 042807.
}
\seealso{\code{\link{compare_layered_models}}, \code{\link{model_evidence}},
\code{\link{graph_set_edge_property}}, \code{demo("layered_sbm")}}
\examples{
set.seed(42)
g <- graph_create(10, directed = FALSE)
edges <- rbind(t(utils::combn(1:5, 2)), t(utils::combn(6:10, 2)),
               c(1, 6), c(2, 7), c(3, 8))
lay <- c(rep("within", 20), rep("across", 3))
graph_add_edges(g, edges)
st <- minimize_layered_blockmodel_dl(g, layers = lay, independent = TRUE)
st
nst <- minimize_nested_layered_blockmodel_dl(g, layers = lay,
                                             independent = TRUE)
nst
}
