\name{bipartite_block_state}
\alias{fungus_tree_network}
\alias{bipartite_graph}
\alias{is_bipartite}
\alias{as_incidence_matrix}
\alias{bipartite_block_state}
\alias{bipartite_description_length}
\alias{minimize_bipartite_blockmodel_dl}
\alias{nested_bipartite_block_state}
\alias{nested_bipartite_description_length}
\alias{minimize_nested_bipartite_blockmodel_dl}
\alias{bipartite_membership}
\alias{bipartite_edge_counts}
\alias{bipartite_block_matrix}
\alias{compare_bipartite_models}
\alias{print.rgraphtool_bipartite_block_state}
\alias{print.rgraphtool_nested_bipartite_block_state}
\alias{print.rgraphtool_bipartite_model_comparison}
\title{Native bipartite network representation and Stochastic Block Model inference}
\usage{
fungus_tree_network()
bipartite_graph(x, type = NULL, directed = FALSE, dim_labels = NULL,
  weight_mode = c("auto", "binary", "multigraph", "weight"),
  weight_name = "weight", covariates = NULL)
is_bipartite(graph, type = NULL)
as_incidence_matrix(graph, type = NULL, weight = NULL)
bipartite_block_state(graph, type = NULL, blocks = NULL,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  edge_prior = c("bipartite", "unipartite"), weights = NULL, types = NULL,
  params = NULL, microcanonical = FALSE, dim_labels = NULL)
bipartite_description_length(state)
minimize_bipartite_blockmodel_dl(graph, type = NULL,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  edge_prior = c("bipartite", "unipartite"), weights = NULL, types = NULL,
  params = NULL, microcanonical = FALSE, B_row_range = NULL, B_col_range = NULL,
  nruns = 1L, refine = TRUE, control = list(), init_blocks = NULL,
  dim_labels = NULL,
  weight_mode = c("auto", "binary", "multigraph", "weight"))
nested_bipartite_block_state(graph, type = NULL, levels = NULL,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  edge_prior = c("bipartite", "unipartite"), weights = NULL, types = NULL,
  params = NULL, microcanonical = FALSE, dim_labels = NULL)
nested_bipartite_description_length(state)
minimize_nested_bipartite_blockmodel_dl(graph, type = NULL,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  edge_prior = c("bipartite", "unipartite"), weights = NULL, types = NULL,
  params = NULL, microcanonical = FALSE, nruns = 1L, control = list(),
  dim_labels = NULL,
  weight_mode = c("auto", "binary", "multigraph", "weight"))
bipartite_membership(state, level = 1L)
bipartite_edge_counts(state, level = 1L)
bipartite_block_matrix(state, level = 1L,
  normalize = c("density", "counts", "mean_weight"))
compare_bipartite_models(graph, type = NULL, weights = NULL, types = NULL,
  params = NULL, microcanonical = FALSE,
  edge_prior = c("bipartite", "unipartite"), nested = TRUE, nruns = 3L,
  dim_labels = NULL)
\method{print}{rgraphtool_bipartite_block_state}(x, ...)
\method{print}{rgraphtool_nested_bipartite_block_state}(x, ...)
\method{print}{rgraphtool_bipartite_model_comparison}(x, ...)
}
\arguments{
  \item{x}{For \code{bipartite_graph()}, a rectangular numeric/logical incidence
  matrix (or sparse \code{Matrix}), an \code{rgraphtool_graph}, or an
  \code{igraph} object with a \code{type} vertex attribute. For \code{print()},
  a bipartite state or comparison object.}
  \item{graph}{An \code{rgraphtool_graph}, rectangular incidence matrix, or
  \code{igraph} object.}
  \item{type}{Optional length-\eqn{N} vertex type vector (logical, 0/1, 1/2,
  character, or factor with at most two distinct values). If \code{NULL}, reads
  vertex property \code{"type"} from \code{graph}.}
  \item{directed}{Logical; whether the constructed bipartite graph is directed
  (default \code{FALSE}).}
  \item{dim_labels}{Optional length-2 character vector naming the row and column
  vertex modes (e.g., \code{c("plant", "pollinator")}).}
  \item{weight_mode}{How nonzero entries of an incidence matrix are encoded:
  \code{"binary"} (presence/absence), \code{"multigraph"} (non-negative integer
  parallel edges), \code{"weight"} (single edges with numeric edge property
  \code{weight_name}), or \code{"auto"}.}
  \item{weight_name}{Edge property name used when \code{weight_mode = "weight"}.}
  \item{covariates}{Optional named list of matrices with \code{dim(x)} to attach
  as additional edge properties.}
  \item{weight}{Optional edge property name to aggregate in
  \code{as_incidence_matrix()}.}
  \item{blocks}{Optional length-\eqn{N} positive integer block assignment
  respecting bipartite vertex types.}
  \item{levels}{Optional list of per-level integer block vectors for
  \code{nested_bipartite_block_state()}.}
  \item{degree_corrected}{Logical; whether to fit the degree-corrected (DC) or
  non-degree-corrected (NDC) microcanonical bipartite SBM.}
  \item{degree_prior}{Degree sequence hyperprior (\code{"distributed"} or
  \code{"uniform"}).}
  \item{edge_prior}{Block edge-count matrix prior: \code{"bipartite"} (default,
  exact Yen & Larremore 2020 \eqn{\ln\binom{B_1 B_2 + E - 1}{E}} prior
  exploiting bipartite structural zeros) or \code{"unipartite"} (matches
  \code{graph-tool} 3.9 \code{BlockState(..., pclabel=type)}).}
  \item{weights}{Optional edge weight property name(s), vector, or matrix for
  weighted bipartite SBMs.}
  \item{types}{Optional covariate distribution type(s) (\code{"real-exponential"},
  \code{"real-normal"}, \code{"discrete-geometric"}, \code{"discrete-poisson"},
  \code{"discrete-binomial"}).}
  \item{params}{Optional list of conjugate hyperparameter lists.}
  \item{microcanonical}{Logical; whether to use microcanonical covariate priors.}
  \item{B_row_range, B_col_range}{Optional length-2 integer bounds
  \code{c(min, max)} on the number of row and column blocks.}
  \item{nruns}{Number of independent restarts; keeps the lowest description
  length.}
  \item{refine}{Logical; whether to perform zero-temperature refinement sweeps.}
  \item{control}{Optional list of search controls (\code{r} shrink ratio,
  \code{n_cand} merge candidates per block).}
  \item{init_blocks}{Optional warm-start block assignment vector.}
  \item{state}{A bipartite or nested bipartite block state.}
  \item{level}{Hierarchy level (1-based) for nested bipartite states.}
  \item{normalize}{Summary returned by \code{bipartite_block_matrix()}:
  \code{"density"} (\eqn{e_{rs} / (n_r n_s)}), \code{"counts"} (\eqn{e_{rs}}),
  or \code{"mean_weight"}.}
  \item{nested}{Logical; whether \code{compare_bipartite_models()} also fits
  nested bipartite models.}
  \item{...}{Unused; included for S3 \code{print} compatibility.}
}
\description{
Construct, validate, infer, and summarize flat, hierarchical (nested), and
weighted microcanonical Bipartite Stochastic Block Models (biSBM) in native
multi-threaded C++23.
}
\details{
The native C++ engine implements the exact microcanonical bipartite description
length of Yen & Larremore (2020, \emph{Phys. Rev. E} 102, 032309) with
type-constrained partitions, exact per-type partition priors, and the bipartite
edge-count prior \eqn{\ln\binom{B_1 B_2 + E - 1}{E}} (or
\eqn{\ln\binom{2 B_1 B_2 + E - 1}{E}} for directed bipartite networks). Setting
\code{edge_prior = "unipartite"} reproduces \code{graph-tool} 3.9's
\code{BlockState(..., pclabel=type)} to machine precision.

Flat bipartite states also integrate with \code{\link{mcmc_sweep}()},
\code{\link{mcmc_equilibrate}()}, and \code{\link{sample_posterior}()} for
exact finite-temperature Bayesian posterior sampling.
}
\value{
\code{fungus_tree_network()} returns a list containing the antagonistic
tree--fungus interaction dataset of Vacher, Piou & Desprez-Loustau (2008,
\emph{PLoS ONE} 3(3): e1740; packaged in \pkg{sbm} as
\code{fungusTreeNetwork}): \code{fungus_tree} (\eqn{154 \times 51} binary
incidence matrix), \code{tree_tree} (\eqn{51 \times 51} shared-fungus count
matrix), \code{covar_tree} (genetic, taxonomic, and geographic distance
matrices), \code{fungus_names}, \code{tree_names}, and \code{graph} (the
bipartite \code{rgraphtool_graph}).
\code{bipartite_graph()} returns an \code{rgraphtool_graph} with vertex property
\code{"type"}. \code{is_bipartite()} returns a logical scalar.
\code{as_incidence_matrix()}, \code{bipartite_edge_counts()}, and
\code{bipartite_block_matrix()} return rectangular numeric matrices.
\code{bipartite_block_state()} and \code{minimize_bipartite_blockmodel_dl()}
return an \code{rgraphtool_bipartite_block_state}; their nested counterparts
return an \code{rgraphtool_nested_bipartite_block_state}.
\code{bipartite_membership()} returns a named length-2 list of integer block
memberships for the two modes.
}
\examples{
set.seed(1)
M <- matrix(0L, 12, 12)
M[1:6, 1:6] <- rbinom(36, 1, 0.85)
M[7:12, 7:12] <- rbinom(36, 1, 0.85)
M[1:6, 7:12] <- rbinom(36, 1, 0.05)
M[7:12, 1:6] <- rbinom(36, 1, 0.05)

fit <- minimize_bipartite_blockmodel_dl(M, dim_labels = c("plant", "pollinator"))
fit
bipartite_membership(fit)
bipartite_block_matrix(fit, normalize = "density")
}
