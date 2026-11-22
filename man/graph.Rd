\name{graph_create}
\alias{graph_create}
\alias{graph_add_vertices}
\alias{graph_add_edges}
\alias{graph_num_vertices}
\alias{graph_num_edges}
\alias{graph_edges}
\alias{graph_degree}
\title{Create and inspect a graph}
\usage{
graph_create(n = 0L, directed = TRUE)
graph_add_vertices(graph, n = 1L)
graph_add_edges(graph, edges)
graph_num_vertices(graph)
graph_num_edges(graph)
graph_edges(graph)
graph_degree(graph, mode = c("total", "in", "out"))
}
\arguments{
  \item{n}{Non-negative number of vertices.}
  \item{directed}{Whether the graph is directed.}
  \item{graph}{A graph created by \code{graph_create()}.}
  \item{edges}{Two-column numeric matrix of 1-based source and target vertex IDs.}
  \item{mode}{Degree type. Undirected graphs support only total degree.}
}
\description{Create and modify a native graph-tool adjacency structure and inspect its contents.}
\details{
Graphs are mutable and aliases share the same underlying graph. External pointers
cannot be restored from R serialization. Parallel edges and self-loops are
allowed; a self-loop contributes two to total degree. Edge endpoints are
returned in their original orientation even for undirected graphs.
}
\value{
\code{graph_create()} returns a graph. Modification functions invisibly return
the graph. Counts return numbers, \code{graph_edges()} returns a two-column
integer matrix, and \code{graph_degree()} returns a numeric vector.
}
\examples{
g <- graph_create(3)
graph_add_edges(g, rbind(c(1, 2), c(2, 3)))
graph_degree(g, "out")
}
