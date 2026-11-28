\name{description_length}
\alias{description_length}
\alias{state_blocks}
\alias{is_nested_state}
\alias{is_block_state}
\title{Common accessors for block states}
\usage{
description_length(state, ...)
state_blocks(state, level = 0L)
is_nested_state(state)
is_block_state(x)
}
\arguments{
  \item{state}{A block state of any family, e.g. from \code{\link{block_state}},
  \code{\link{nested_block_state}}, \code{\link{weighted_block_state}},
  \code{\link{layered_block_state}}, \code{\link{pp_block_state}} or a
  minimizer.}
  \item{...}{Passed to the family's own description-length function (for
  example \code{degree_prior} or \code{lmemb}).}
  \item{level}{Hierarchy level (0 = finest). Flat states have only level 0;
  dynamics states use the hierarchy of their SBM prior.}
  \item{x}{Any object.}
}
\description{
\code{description_length()} is an S3 generic returning a state's description
length in nats; each method calls the family's own function
(\code{block_description_length()}, \code{nested_description_length()},
\code{weighted_description_length()}, \code{layered_description_length()},
\code{pp_description_length()}, \code{dynamics_description_length()}, and so on).
\code{state_blocks()} returns the vertex partition at a hierarchy level for flat
and nested states (and the SBM prior of dynamics states).
\code{is_nested_state()} tells flat and nested states apart, and
\code{is_block_state()} is \code{TRUE} for a state of any family, so code
handling several families needs no class checks of its own.
}
\value{
\code{description_length()}: a numeric scalar. \code{state_blocks()}: an
integer vector with one group label per vertex. \code{is_nested_state()},
\code{is_block_state()}: a logical scalar.
}
\seealso{\code{\link{block_membership}}, \code{\link{nested_project}}}
\examples{
g <- graph_create(6, directed = FALSE)
graph_add_edges(g, rbind(c(1, 2), c(2, 3), c(1, 3), c(4, 5), c(5, 6), c(4, 6)))
st <- block_state(g, c(1, 1, 1, 2, 2, 2))
description_length(st)
state_blocks(st)
nst <- nested_block_state(g, list(c(1, 1, 1, 2, 2, 2), c(1, 1)))
description_length(nst)
state_blocks(nst, level = 1)
is_nested_state(nst)
}
