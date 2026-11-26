\name{expression_utils}
\alias{project_fixed_effects}
\alias{filter_expression_genes}
\alias{module_eigengenes}
\alias{shrinkage_covariance}
\title{Gene-Expression Fixed-Effects Projection, Filtering, Shrinkage Covariance, and Module Eigengenes}
\description{
Utilities for conditioning gene-expression matrices on fixed-effect designs via
exact QR/REML orthogonal projection, pre-filtering informative genes, computing
analytical shrinkage covariances, and computing sign-aligned module eigengenes.
}
\usage{
project_fixed_effects(
  expr,
  design,
  col_data = NULL,
  transform = c("zscore", "nonparanormal", "none")
)

filter_expression_genes(
  expr,
  col_data = NULL,
  nuisance = NULL,
  signal = NULL,
  min_var = 0,
  max_padj = 1,
  top_n = NULL
)

module_eigengenes(expr, blocks, level = 0L, min_size = 2L)

shrinkage_covariance(expr, target = c("diagonal", "identity"))
}
\arguments{
\item{expr}{Numeric matrix of genes (rows) by samples (columns).}
\item{design}{A formula evaluated on \code{col_data} or an \code{M x p} numeric design matrix.}
\item{col_data}{Optional \code{data.frame} of sample metadata.}
\item{transform}{One of \code{"zscore"}, \code{"nonparanormal"}, or \code{"none"}.}
\item{nuisance}{Optional nuisance formula or matrix.}
\item{signal}{Optional signal formula or matrix tested beyond \code{nuisance}.}
\item{min_var}{Minimum residual variance after removing \code{nuisance}.}
\item{max_padj}{Maximum Benjamini-Hochberg adjusted p-value for \code{signal}.}
\item{top_n}{Optional maximum number of genes to retain.}
\item{blocks}{Integer block vector or \code{rgraphtool} block state.}
\item{level}{Hierarchy level (\code{0L} = finest) for nested states.}
\item{min_size}{Minimum block size for eigengene calculation.}
\item{target}{Shrinkage target (\code{"diagonal"} or \code{"identity"}).}
}
\value{
\code{project_fixed_effects()} returns an object of class
\code{"rgraphtool_expression_projection"}; \code{filter_expression_genes()}
returns a list with \code{expr}, \code{stats}, and \code{keep};
\code{module_eigengenes()} returns a list with \code{eigengenes}, \code{pve},
and \code{block_sizes}; \code{shrinkage_covariance()} returns a list with
\code{cov}, \code{lambda}, and \code{n_obs}.
}
