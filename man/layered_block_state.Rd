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
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"))
layered_description_length(state, lmemb = TRUE, degree_prior = NULL)
minimize_layered_blockmodel_dl(graph, layers, independent = TRUE,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  B_min = 1L, B_max = Inf, refine = TRUE, epsilon = 1e-3)
nested_layered_block_state(graph, layers, independent = TRUE, levels = NULL,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"))
nested_layered_level_entropy(state, degree_prior = NULL)
nested_layered_description_length(state, degree_prior = NULL)
minimize_nested_layered_blockmodel_dl(graph, layers, independent = TRUE,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  B_min = 1L, B_max = Inf, refine = TRUE, epsilon = 1e-3)
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
  (applicable when \code{independent = TRUE}).}
  \item{B_min, B_max}{Minimum and maximum number of occupied level-0 groups.}
  \item{refine}{Whether to run local refinement sweeps during minimization.}
  \item{epsilon}{Convergence tolerance in nats.}
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
