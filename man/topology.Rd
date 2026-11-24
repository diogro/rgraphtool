\name{shortest_distance}
\alias{shortest_distance}
\alias{shortest_path}
\alias{pseudo_diameter}
\alias{graph_diameter}
\alias{label_components}
\alias{label_biconnected_components}
\alias{topological_sort}
\alias{is_dag}
\alias{local_clustering}
\alias{global_clustering}
\alias{assortativity}
\title{Shortest paths, connected components, clustering, and assortativity}
\usage{
shortest_distance(graph, source = NULL, target = NULL, weights = NULL,
                  max_dist = Inf)
shortest_path(graph, source, target, weights = NULL)
pseudo_diameter(graph, source = 1L, weights = NULL)
graph_diameter(graph, weights = NULL)
label_components(graph)
label_biconnected_components(graph)
topological_sort(graph)
is_dag(graph)
local_clustering(graph)
global_clustering(graph)
assortativity(graph, values = "total", categorical = FALSE, weights = NULL)
}
\arguments{
  \item{graph}{An \code{rgraphtool_graph} object.}
  \item{source}{Optional 1-based source vertex ID(s). When \code{NULL} in
  \code{shortest_distance()}, computes all-pairs distances in parallel.}
  \item{target}{Optional 1-based target vertex ID(s).}
  \item{weights}{Optional edge property name or numeric vector of edge lengths /
  weights.}
  \item{max_dist}{Maximum distance cutoff for \code{shortest_distance()}.}
  \item{values}{For \code{assortativity()}, \code{"total"}, \code{"in"}, or
  \code{"out"} for degree assortativity, a vertex property name, or a length-\eqn{N}
  vector.}
  \item{categorical}{Logical; if \code{TRUE}, computes discrete Newman
  assortativity rather than scalar Pearson correlation.}
}
\description{
Native C++23 graph topology and descriptive network algorithms: unweighted
(BFS), non-negative weighted (Dijkstra), and negative-weight (Bellman--Ford)
shortest distances and path reconstruction, exact and pseudo-diameter,
connected and strongly connected components (with attractor detection),
Hopcroft--Tarjan biconnected components (articulation points and bridges),
topological sorting and DAG testing, local and global clustering coefficients
(with jackknife standard error), and scalar/categorical assortativity.
}
\value{
\code{shortest_distance()} returns a numeric vector (for a single
\code{source}) or matrix. \code{shortest_path()} returns
\code{list(vertices, edges, distance)}. \code{pseudo_diameter()} returns
\code{list(diameter, endpoints)} and \code{graph_diameter()} returns a scalar.
\code{label_components()} returns \code{list(membership, sizes, num_components,
attractors)}. \code{label_biconnected_components()} returns
\code{list(edge_component, num_components, articulation_points, bridges)}.
\code{topological_sort()} returns an integer vertex ordering and
\code{is_dag()} returns a logical scalar. \code{local_clustering()} returns a
length-\eqn{N} vector, \code{global_clustering()} returns a named length-2
vector \code{c(clustering, std_err)}, and \code{assortativity()} returns a
numeric correlation coefficient.
}
\examples{
g <- graph_create(6L, directed = FALSE)
graph_add_edges(g, rbind(c(1, 2), c(2, 3), c(3, 1), c(3, 4), c(4, 5), c(5, 6), c(6, 4)))
shortest_path(g, 1L, 6L)
label_biconnected_components(g)
global_clustering(g)
assortativity(g, "total")
}
