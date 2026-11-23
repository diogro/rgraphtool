\name{predict_edges}
\alias{predict_edges}
\alias{evaluate_edge_prediction}
\alias{print.rgraphtool_edge_prediction}
\alias{print.rgraphtool_edge_evaluation}
\title{Bayesian missing and spurious edge prediction and evaluation}
\usage{
predict_edges(x, candidates, nsamples = 100L, burnin = 50L, thin = 1L,
  nested = TRUE, degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"),
  self_loops = FALSE, epsilon = 1e-8)

evaluate_edge_prediction(scores, labels, threshold = 0.5)
}
\arguments{
  \item{x}{A \code{rgraphtool_graph}, a fitted \code{rgraphtool_block_state} or
  \code{rgraphtool_nested_block_state}, or a
  \code{rgraphtool_uncertain_block_state}.}
  \item{candidates}{Two-column integer matrix or data.frame of 1-based candidate
  vertex pairs \code{(u, v)} to score.}
  \item{nsamples, burnin, thin}{Posterior MCMC sample count, burn-in sweeps, and
  thinning interval.}
  \item{nested, degree_corrected, degree_prior}{SBM prior options used when
  \code{x} is a raw graph.}
  \item{self_loops}{Whether self-loops are valid candidate dyads.}
  \item{epsilon}{Multiplicity summation convergence tolerance.}
  \item{scores}{Numeric vector of predicted edge scores/probabilities, or a
  \code{rgraphtool_edge_prediction} object.}
  \item{labels}{Binary ground-truth labels (\code{0}/\code{1} or logical) with
  at least one positive and one negative example.}
  \item{threshold}{Classification threshold for precision, recall, F1, and
  accuracy.}
}
\description{
Computes model-averaged Bayesian posterior edge probabilities for candidate
missing or spurious links using flat/nested SBMs or uncertain-network models,
and evaluates binary link-prediction performance with exact tie-aware ROC-AUC,
PR-AUC, Brier score, and ROC/PR curves.
}
\value{
\code{predict_edges()} returns a \code{rgraphtool_edge_prediction} containing a
\code{predictions} data.frame (\code{from}, \code{to}, \code{prob},
\code{log_prob}, and \code{observed} or \code{conditional_prob}).

\code{evaluate_edge_prediction()} returns a \code{rgraphtool_edge_evaluation}
containing \code{roc_auc}, \code{pr_auc}, \code{precision}, \code{recall},
\code{f1}, \code{accuracy}, \code{brier}, \code{confusion}, and full
\code{roc} and \code{pr} curve data.frames.
}
\references{
Peixoto, T. P. (2018). Reconstructing networks with unknown and heterogeneous
errors. \emph{Physical Review X} 8, 041011. doi:10.1103/PhysRevX.8.041011
}
\seealso{\code{\link{measured_block_state}}, \code{\link{edge_probabilities}},
\code{demo("uncertain_reconstruction")}}
\examples{
g <- graph_create(12, directed = FALSE)
e <- rbind(t(combn(1:6, 2)), t(combn(7:12, 2)))
graph_add_edges(g, e[-1, , drop = FALSE])
cands <- rbind(e[1, ], c(1L, 12L))
set.seed(1)
pred <- predict_edges(g, cands, nsamples = 20L, burnin = 10L)
evaluate_edge_prediction(pred, c(1L, 0L))
}
