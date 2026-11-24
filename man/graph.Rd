\name{graph_create}
\alias{graph_create}
\alias{graph_add_vertices}
\alias{graph_add_edges}
\alias{graph_num_vertices}
\alias{graph_num_edges}
\alias{graph_edges}
\alias{graph_degree}
\alias{graph_set_edge_property}
\alias{graph_edge_property}
\alias{graph_edge_properties}
\alias{graph_set_vertex_property}
\alias{graph_vertex_property}
\alias{graph_vertex_properties}
\alias{graph_set_graph_property}
\alias{graph_graph_property}
\alias{graph_graph_properties}
\alias{print.rgraphtool_graph}
\alias{transform_weights}
\title{Create, modify, and inspect a graph and its property maps}
\usage{
graph_create(n = 0L, directed = TRUE)
graph_add_vertices(graph, n = 1L)
graph_add_edges(graph, edges, properties = NULL)
graph_num_vertices(graph)
graph_num_edges(graph)
graph_edges(graph)
graph_degree(graph, mode = c("total", "in", "out"))
graph_set_edge_property(graph, name, values)
graph_edge_property(graph, name)
graph_edge_properties(graph)
graph_set_vertex_property(graph, name, values)
graph_vertex_property(graph, name)
graph_vertex_properties(graph)
graph_set_graph_property(graph, name, value)
graph_graph_property(graph, name)
graph_graph_properties(graph)
\method{print}{rgraphtool_graph}(x, max_edges = 20L, max_ascii = 12L,
  ascii = NULL, ...)
transform_weights(x, transform = c("identity", "log", "log1p", "sqrt",
  "sign", "abs", "rank", "zscore"), shift = 0)
}
\arguments{
  \item{n}{Non-negative number of vertices.}
  \item{directed}{Whether the graph is directed.}
  \item{graph}{A graph created by \code{graph_create()}.}
  \item{edges}{Two-column numeric matrix of 1-based source and target vertex IDs.}
  \item{properties}{Optional named list or data frame of edge vectors
  (each of length \code{nrow(edges)}) to attach as edge property maps aligned
  with the newly added edges.}
  \item{mode}{Degree type. Undirected graphs support only total degree.}
  \item{name}{Non-empty character string naming a graph, vertex, or edge property map.}
  \item{values}{Vector, matrix, or list of property values aligned with
  \code{graph_edges(graph)} (for edge properties) or \code{1:graph_num_vertices(graph)}
  (for vertex properties).}
  \item{value}{Single value (scalar, vector, or raw object) for a graph-level property.}
  \item{x}{For \code{print()}, an \code{rgraphtool_graph} object; for
  \code{transform_weights()}, a finite numeric vector of edge weights to transform.}
  \item{max_edges}{Maximum number of edges to display in the edge preview.}
  \item{max_ascii}{Maximum vertex count for automatically displaying the ASCII
  adjacency matrix.}
  \item{ascii}{Optional logical overriding whether to display the ASCII
  adjacency matrix.}
  \item{...}{Unused; included for S3 \code{print} compatibility.}
  \item{transform}{Weight transformation to apply. Invertible continuous
  transformations (\code{"identity"}, \code{"log"}, \code{"log1p"},
  \code{"sqrt"}, \code{"zscore"}) attach the exact log-Jacobian
  \eqn{\sum_e \ln |f'(x_e)|} as attribute \code{"log_jacobian"} for
  change-of-variables model comparison.}
  \item{shift}{Optional finite numeric constant added to \code{x} before
  applying \code{transform}.}
}
\description{
Create and modify a native graph-tool adjacency structure, inspect its vertices
and edges, attach numeric edge property maps, and transform edge weights with
exact change-of-variables log-Jacobians.
}
\details{
Graphs are mutable and aliases share the same underlying graph. External pointers
cannot be restored from R serialization. Parallel edges and self-loops are
allowed; a self-loop contributes two to total degree. Edge endpoints are
returned in their original orientation even for undirected graphs.

Edge properties passed through \code{graph_add_edges(properties = ...)} are
stored by internal edge insertion index and returned by
\code{graph_edge_property()} in the exact row order of \code{graph_edges(graph)},
even when edges are inserted out of vertex order.
}
\value{
\code{graph_create()} returns a graph. Modification functions invisibly return
the graph. Counts return numbers, \code{graph_edges()} returns a two-column
integer matrix, \code{graph_degree()} and \code{graph_edge_property()} return
numeric vectors, \code{graph_edge_properties()} returns a character vector of
property names, and \code{transform_weights()} returns a numeric vector with
attributes \code{"transform"} and \code{"log_jacobian"}.
}
\examples{
g <- graph_create(3)
graph_add_edges(g, rbind(c(1, 2), c(2, 3)),
                properties = list(weight = c(1.5, 4.2)))
graph_degree(g, "out")
w <- graph_edge_property(g, "weight")
transform_weights(w, "log")
}
