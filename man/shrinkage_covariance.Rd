\name{shrinkage_covariance}
\alias{shrinkage_covariance}
\title{Analytical Ledoit-Wolf / Schafer-Strimmer Shrinkage Covariance}
\description{
Computes an optimal linear shrinkage covariance matrix
\eqn{(1 - \lambda) S + \lambda T} for high-dimensional data
(\eqn{N \gg M}) with zero external package dependencies, suitable for passing
directly to \code{\link{pseudo_normal_block_state}(cov = ..., n_eff = ...)}.
}
\usage{
shrinkage_covariance(x, target = c("diagonal", "identity"))
}
\arguments{
\item{x}{Numeric \code{N x M} matrix of variables (rows) by observations
(columns), with \code{M >= 3}.}
\item{target}{Shrinkage target: \code{"diagonal"} (unequal variable variances,
zero covariances) or \code{"identity"} (mean-variance scaled identity).}
}
\value{
A list with \code{cov} (\code{N x N} symmetric positive-definite matrix),
\code{lambda} (shrinkage intensity in \code{[0, 1]}), and \code{n_obs} (\code{M}).
}
\seealso{
\code{\link{pseudo_normal_block_state}}
}
\examples{
set.seed(1)
x <- matrix(rnorm(20 * 10), nrow = 20, ncol = 10)
sh <- shrinkage_covariance(x, target = "diagonal")
sh$lambda
dim(sh$cov)
}
