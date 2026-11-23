\name{rgraphtool_threads}
\alias{rgraphtool_max_threads}
\alias{rgraphtool_num_threads}
\alias{rgraphtool_set_num_threads}
\title{Configure parallel worker threads for large-network inference}
\usage{
rgraphtool_max_threads()
rgraphtool_num_threads()
rgraphtool_set_num_threads(nthreads)
}
\arguments{
  \item{nthreads}{Positive integer number of worker threads (capped at 256), or
  \code{"auto"} to use \code{rgraphtool_max_threads()}.}
}
\description{
Query and configure the portable C++23 thread pool used to parallelize
large-network stochastic block model minimization in
\code{\link{minimize_blockmodel_dl}}, \code{\link{refine_blockmodel_dl}}, and
\code{\link{minimize_nested_blockmodel_dl}}.
}
\details{
\pkg{rgraphtool} uses a built-in C++23 \code{std::thread} pool rather than
OpenMP, so multithreading works identically on macOS (including Apple clang
without \code{libomp}), Linux, and Windows with no external threading library.
Worker threads never call the R API or R's RNG; per-thread deterministic
\code{xoshiro256**} generators are seeded from R's RNG on the main thread before
entering parallel regions, so results remain reproducible under \code{set.seed()}
for a fixed thread count.

By default \pkg{rgraphtool} runs on \code{1L} thread so background tasks and
shared systems are never oversubscribed unexpectedly. Opt in globally via
\code{rgraphtool_set_num_threads(k)} or \code{options(rgraphtool.nthreads = k)},
or per call via the \code{nthreads} argument.

In practice, 4--16 threads typically provide the best speedups (roughly 3x--4.3x
on networks with \eqn{N = 10^4}--\eqn{10^5} vertices), after which sequential
priority-queue merge commits limit further scaling; 32--64 threads still give
incremental gains on \eqn{N \ge 10^5} networks. To reduce peak RAM on very large
networks, lower \code{options(rgraphtool.q_cache_max = 1000)} (the default
\code{5000} uses ~100 MB for the exact partition-count table) and pair
multithreading with strictly agglomerative search
(\code{control = list(bisection = FALSE, r = 0.5)}) or hierarchy-shape
initialization (\code{simple_init = TRUE}).
}
\value{
\code{rgraphtool_max_threads()} returns the number of hardware threads reported
by \code{std::thread::hardware_concurrency()} (at least \code{1L}).

\code{rgraphtool_num_threads()} returns the active default thread count.

\code{rgraphtool_set_num_threads()} invisibly returns the configured integer
thread count.
}
\seealso{\code{\link{minimize_blockmodel_dl}},
\code{\link{minimize_nested_blockmodel_dl}}, \code{demo("large_networks")}}
\examples{
rgraphtool_max_threads()
old <- rgraphtool_num_threads()
rgraphtool_set_num_threads(min(2L, rgraphtool_max_threads()))
rgraphtool_num_threads()
rgraphtool_set_num_threads(old)
}
