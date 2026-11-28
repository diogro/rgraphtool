\name{shrink_pcor}
\alias{shrink_pcor}
\alias{shrink_pcor_lambda}
\title{Fast Shrinkage Partial Correlations}
\description{
Schafer--Strimmer shrinkage partial correlations: the sample correlation
matrix is shrunk towards the identity,
\eqn{R^* = (1 - \lambda) R + \lambda I}, and inverted. The inverse is
computed through the Woodbury identity in the \eqn{M \times M} observation
space, so for \eqn{N \gg M} no \eqn{N \times N} matrix is formed unless
requested: partial correlations can be returned for selected pairs, or as
the \code{top_k} strongest pairs, streaming over blocks of variables.
}
\usage{
shrink_pcor(x, lambda = NULL, center = TRUE, pairs = NULL, top_k = NULL,
            matrix = NULL)

shrink_pcor_lambda(x, center = TRUE)
}
\arguments{
\item{x}{Numeric \code{N x M} matrix of variables (rows) by observations
(columns), \code{M >= 3}. Rows must have non-zero variance.}
\item{lambda}{Shrinkage intensity in \code{[0, 1]}, or \code{NULL} (default)
to estimate the Schafer--Strimmer optimal intensity from \code{x}. Fixing
\code{lambda} keeps estimates from different data sets (for example two
conditions of different size) on the same shrinkage scale.}
\item{center}{Logical; centre each variable (default \code{TRUE}). With
\code{FALSE} the observations are taken to have known zero mean (for
example regression residuals), correlations use \eqn{\sum_k x_{ik} x_{jk}}
without centring, and the intensity estimate uses the matching variance
formula.}
\item{pairs}{Optional two-column integer matrix of 1-based variable pairs at
which to return partial correlations.}
\item{top_k}{Optional number of strongest pairs (by absolute partial
correlation, \code{i < j}) to return.}
\item{matrix}{Logical; return the full \code{N x N} matrix. Defaults to
\code{TRUE} only when neither \code{pairs} nor \code{top_k} is given.}
}
\value{
\code{shrink_pcor_lambda()} returns the estimated intensity only (cost
\eqn{O(N M^2)}). \code{shrink_pcor()} returns a list with \code{lambda} (the intensity used), \code{n_obs} (\code{M}), and
as requested \code{matrix} (\code{N x N}, zero diagonal), \code{pairs}
(numeric vector aligned with \code{pairs}), and \code{top} (data frame with
\code{from}, \code{to}, \code{pcor}, sorted by decreasing
\code{abs(pcor)}).
}
\details{
With \eqn{A} the standardized data scaled so that \eqn{R = A A^\top},
\eqn{c = \lambda / (1 - \lambda)} and the Cholesky factorization
\eqn{c I_M + A^\top A = L^\top L}, the partial correlations are
\eqn{\rho_{ij} = v_i^\top v_j} with \eqn{u_i} the rows of \eqn{A L^{-1}} and
\eqn{v_i = u_i / \sqrt{1 - |u_i|^2}}. Time is \eqn{O(N M^2)} for the
intensity and factorization plus \eqn{O(N^2 M)} for \code{top_k} or the full
matrix; memory is \eqn{O(N M + M^2)} plus the requested output. With
\code{lambda = 0} the correlation matrix is inverted directly and must be
non-singular. Results match \code{corpcor::pcor.shrink()} for
\code{center = TRUE}.
}
\references{
Schafer, J. and Strimmer, K. (2005). A shrinkage approach to large-scale
covariance matrix estimation and implications for functional genomics.
\emph{Statistical Applications in Genetics and Molecular Biology} 4(1):32.
}
\seealso{
\code{\link{shrinkage_covariance}}
}
\examples{
set.seed(1)
x <- matrix(rnorm(50 * 20), nrow = 50, ncol = 20)
pc <- shrink_pcor(x)
pc$lambda
top <- shrink_pcor(x, top_k = 5)$top
top
}
