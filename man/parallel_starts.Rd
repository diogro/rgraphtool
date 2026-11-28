\name{parallel_starts}
\alias{parallel_starts}
\title{Run Independent Starts in Parallel Processes}
\description{
Evaluates \code{fun(1)}, ..., \code{fun(n)} -- typically independent
minimizations from different random starts -- in up to \code{cores} forked
worker processes, and returns the list of results. This is the engine behind
the \code{cores} argument of the minimizers and model comparisons that take
\code{nruns}.
}
\usage{
parallel_starts(n, fun, cores = NULL, graph = NULL)
}
\arguments{
\item{n}{Number of starts.}
\item{fun}{Function of the run index returning one result (for example a
fitted state).}
\item{cores}{Number of worker processes; \code{NULL} uses
\code{getOption("rgraphtool.cores", 1L)}. Forked workers are unavailable on
Windows, where starts run sequentially with a warning.}
\item{graph}{Optional \code{rgraphtool_graph}. External pointers do not
survive the return from a worker, so every graph object inside a result is
replaced by \code{graph}.}
}
\value{A list of the \code{n} results, in run order.}
\details{
With \code{n > 1}, a seed is drawn for every run before any run starts and
each run sets its own seed, so results are identical for any value of
\code{cores} (and the caller's random stream advances the same way). With
\code{n = 1} and \code{cores = 1} the function simply calls \code{fun(1)}.
Workers run with one rgraphtool thread (the parent's thread pool does not
exist in a forked child); total CPU use is therefore \code{cores} processes.
Set \code{OPENBLAS_NUM_THREADS=1} when using a multithreaded BLAS.
}
\examples{
g <- graph_create(30, directed = FALSE)
set.seed(1)
graph_add_edges(g, cbind(sample(30, 90, TRUE), sample(30, 90, TRUE)))
fits <- parallel_starts(3, function(run) minimize_blockmodel_dl(g), graph = g)
sapply(fits, block_description_length)
}
