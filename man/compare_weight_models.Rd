\name{compare_weight_models}
\alias{compare_weight_models}
\alias{print.rgraphtool_weight_model_comparison}
\title{Compare weighted SBM covariate distributions and transformations}
\usage{
compare_weight_models(graph, weights, candidates = NULL, nested = TRUE,
  degree_corrected = TRUE, degree_prior = c("distributed", "uniform"),
  method = c("map", "posterior"),
  evidence_method = c("variational", "importance"), nruns = 1L, cores = NULL,
  nsamples = 200L, thin = 5L, burnin = NULL, mle = TRUE, B_min = 1L,
  B_max = Inf, refine = TRUE)
}
\arguments{
  \item{graph}{A graph created by \code{\link{graph_create}}.}
  \item{weights}{Single edge weight covariate: an edge property name or numeric
  vector of length \code{graph_num_edges(graph)}.}
  \item{candidates}{Candidate weight models to compare:
  \itemize{
    \item \code{NULL} (default): automatically selects discrete count models
    (\code{"discrete-geometric"}, \code{"discrete-poisson"},
    \code{"discrete-binomial"}) for non-negative integer weights, continuous
    positive models (\code{"real-exponential"}, \code{"real-normal"},
    \code{"log-normal"}) for strictly positive real weights, or
    microcanonical vs Bayesian \code{"real-normal"} for signed real weights.
    \item A character vector of shorthand names from
    \code{"real-exponential"}, \code{"real-normal"}, \code{"log-normal"},
    \code{"discrete-geometric"}, \code{"discrete-poisson"},
    \code{"discrete-binomial"}.
    \item A named list where each entry is a list with elements \code{type},
    and optional \code{transform} (default \code{"identity"}), \code{shift}
    (default \code{0}), and \code{params} (default \code{"microcanonical"}).
  }}
  \item{nested}{Whether to fit nested (\code{TRUE}, default) or flat
  (\code{FALSE}) weighted SBMs.}
  \item{degree_corrected}{Whether to use degree correction.}
  \item{degree_prior}{Degree sequence hyperprior (\code{"distributed"} or
  \code{"uniform"}).}
  \item{method}{Rank candidate models by \code{"map"} description length or
  full Bayesian \code{"posterior"} marginal evidence.}
  \item{evidence_method}{Evidence estimator forwarded to
  \code{\link{model_evidence}} when \code{method = "posterior"}.}
  \item{nruns}{Number of independent minimization restarts per candidate.}
  \item{cores}{Number of worker processes running the \code{nruns} independent starts in parallel (forked; \code{NULL} uses \code{getOption("rgraphtool.cores", 1L)}). Workers run single-threaded; results do not depend on \code{cores}. Unavailable on Windows.}
  \item{nsamples, thin, burnin, mle}{Posterior sampling arguments used when
  \code{method = "posterior"}.}
  \item{B_min, B_max}{Bounds on the number of level-0 groups.}
  \item{refine}{Whether to run refinement sweeps during minimization.}
}
\description{
Declarative model selection across weighted stochastic block model covariate
distributions and weight transformations (such as \code{log-normal} vs
\code{real-exponential} and \code{real-normal}), automatically applying the
exact change-of-variables Jacobian adjustment so transformed models are
compared on the original weight scale.
}
\details{
For an invertible continuous weight transformation \eqn{y_e = f(w_e)}, the
induced density on the original weights is
\eqn{p_W(w) = p_Y(f(w)) \prod_e |f'(w_e)|}, so the adjusted description length
and negative log-evidence are
\deqn{-\ln p_W(w) = -\ln p_Y(f(w)) - \sum_e \ln |f'(w_e)|.}
For \code{log-normal} (\eqn{f(w) = \ln w}), the Jacobian adjustment added to
the transformed description length is \eqn{+\sum_e \ln w_e}.

Note that discrete probability mass functions and continuous probability
density functions live on different reference measures; by default
(\code{candidates = NULL}), discrete weights are compared against discrete
models and continuous weights against continuous models.
}
\value{
A \code{rgraphtool_weight_model_comparison} list containing:
\itemize{
  \item \code{summary}: ranked data frame with \code{model}, \code{type},
  \code{transform}, \code{B}, \code{raw_dl}, \code{jacobian_adj},
  \code{description_length}, \code{neg_log_evidence},
  \code{posterior_entropy}, \code{delta}, and \code{posterior_weight}.
  \item \code{best_model}: name of the winning candidate.
  \item \code{best_state}: fitted weighted state for the winning candidate.
  \item \code{states}: named list of fitted states ordered from best to worst.
  \item \code{evidences}: named list of \code{rgraphtool_model_evidence} objects
  (when \code{method = "posterior"}).
  \item \code{log_bayes_factors}: pairwise matrix of log Bayes factors.
}
}
\references{
Peixoto, T. P. (2018). Nonparametric weighted stochastic block models.
\emph{Physical Review E} 97, 012306.
}
\seealso{\code{\link{weighted_block_state}}, \code{\link{transform_weights}},
\code{\link{model_evidence}}, \code{demo("weighted_sbm")}}
\examples{
set.seed(42)
g <- graph_create(12, directed = FALSE)
edges <- rbind(t(utils::combn(1:6, 2)), t(utils::combn(7:12, 2)), c(6, 7))
w <- c(exp(stats::rnorm(15, 2, 0.8)), exp(stats::rnorm(15, 0, 0.8)), 1.5)
graph_add_edges(g, edges, properties = list(w = w))
cmp <- compare_weight_models(g, weights = "w", nested = FALSE)
cmp
}
