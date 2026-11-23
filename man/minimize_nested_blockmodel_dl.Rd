\name{minimize_nested_blockmodel_dl}
\alias{minimize_nested_blockmodel_dl}
\alias{nested_block_state}
\alias{nested_description_length}
\alias{nested_level_entropy}
\alias{nested_levels}
\alias{nested_project}
\alias{nested_level_state}
\alias{nested_summary}
\alias{print.rgraphtool_nested_block_state}
\title{Hierarchical (nested) stochastic block models}
\usage{
minimize_nested_blockmodel_dl(graph, degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"), B_min = 1, B_max = Inf,
  refine = TRUE, epsilon = 1e-3, control = list())
nested_block_state(graph, levels, degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"))
nested_description_length(state, degree_prior = NULL)
nested_level_entropy(state, degree_prior = NULL)
nested_levels(state)
nested_project(state, level = 0L)
nested_level_state(state, level = 0L)
nested_summary(state)
}
\arguments{
  \item{graph}{A graph created by \code{graph_create()}, with at least one vertex.}
  \item{degree_corrected}{Whether the lowest level uses the degree-corrected model.}
  \item{degree_prior}{Degree-sequence prior of the lowest level. \code{NULL} uses the
  prior stored in \code{state}.}
  \item{B_min, B_max}{Bounds on the number of groups at level 0.}
  \item{refine}{Whether to run the refinement loop after the level-by-level search.}
  \item{epsilon}{Relative change of the description length that stops refinement.}
  \item{control}{Named list of search settings, as in
  \code{\link{minimize_blockmodel_dl}}. Here \code{c} defaults to 1 and is scaled
  by \eqn{2^l} at level \eqn{l}, as in graph-tool.}
  \item{levels}{A list of label vectors. \code{levels[[1]]} labels each vertex;
  \code{levels[[l]]} labels each group of \code{levels[[l - 1]]}, so it has one
  entry per distinct label of the previous level, in order of first appearance.}
  \item{state}{A nested block state.}
  \item{level}{Hierarchy level, counted from 0 (the vertex partition) as in
  graph-tool.}
}
\description{
The nested SBM places the groups of each level into groups at the level above,
up to a single group. Each upper level is a non-degree-corrected SBM of the
multigraph of group-pair edge counts below it. These functions are the R
counterparts of graph-tool's \code{NestedBlockState} and
\code{minimize_nested_blockmodel_dl()}.
}
\details{
The description length is that of graph-tool 3.9's
\code{NestedBlockState.entropy()} with default arguments. Level 0 is the flat SBM
of \code{\link{block_state}} without the edge-count prior. Level \eqn{l \ge 1}
adds the partition prior of its \eqn{B_{l-1}} nodes and the dense multigraph
likelihood \eqn{\sum_{r \le s} \log \left(\!\binom{n_r n_s}{e_{rs}}\!\right)} of
the group-pair counts. The top level adds the edge-count prior.

\code{minimize_nested_blockmodel_dl()} follows graph-tool's algorithm. It starts
with one group at every level and \eqn{\lceil \log_2 N \rceil + 1} levels. It
then runs the multilevel search of \code{\link{minimize_blockmodel_dl}} on one
level at a time, starting at level 0. Every move accounts for the change in
description length at all higher levels. A result is kept only if it lowers the
total description length. After an improvement the search moves one level up;
otherwise it moves one level down, stopping when it falls below level 0.
Refinement then sweeps all levels until the relative change is at most
\code{epsilon}. The differences from graph-tool listed in
\code{\link{minimize_blockmodel_dl}} apply here too. In addition, if the
single-group start violates \code{B_min} or \code{B_max}, the first level-0
result is accepted even when it does not lower the description length.
graph-tool would keep the single group in that case.

\code{nested_level_entropy()} returns each level's contribution; they sum to
\code{nested_description_length()}. \code{nested_project()} maps vertices to
their groups at a given level. \code{nested_level_state()} returns that
projection as a flat \code{\link{block_state}}. \code{nested_summary()} lists,
for each level up to the first single-node level, the number of nodes \eqn{N},
the number of groups \eqn{B}, and its contribution to the description length.
Printing a state shows the same summary as graph-tool's \code{print_summary()}.

Nested states keep a reference to the graph and fail with an error if the graph
gains edges or vertices afterwards.
}
\value{
\code{minimize_nested_blockmodel_dl()} and \code{nested_block_state()} return a
nested block state: a list with \code{levels} (1-based contiguous labels),
\code{graph}, \code{degree_corrected}, \code{degree_prior} and, for fitted
states, \code{fit}. \code{fit} holds \code{description_length} and
\code{trace}, the total description length after each level search and each
refinement pass.
}
\references{
Peixoto, T. P. (2014). Hierarchical block structures and high-resolution model
selection in large networks. \emph{Physical Review X} 4, 011047.
doi:10.1103/PhysRevX.4.011047
}
\seealso{\code{\link{plot_hierarchy}}, \code{\link{minimize_blockmodel_dl}},
\code{demo("celegans")}}
\examples{
g <- graph_create(12, directed = FALSE)
graph_add_edges(g, do.call(rbind, lapply(list(1:4, 5:8, 9:12),
                                         function(v) t(combn(v, 2)))))
set.seed(1)
state <- minimize_nested_blockmodel_dl(g)
state
nested_summary(state)
nested_project(state, 0)

fixed <- nested_block_state(g, list(rep(1:3, each = 4), c(1, 1, 1)))
nested_description_length(fixed)
}
