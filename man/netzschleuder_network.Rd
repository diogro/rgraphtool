\name{netzschleuder_network}
\alias{netzschleuder_network}
\title{Load a network from the Netzschleuder repository}
\usage{
netzschleuder_network(name, network = name, file = NULL, directed = NULL)
}
\arguments{
  \item{name}{Netzschleuder dataset name, for example \code{"celegansneural"}.}
  \item{network}{Network within the dataset; defaults to \code{name}.}
  \item{file}{Optional path to a previously downloaded \code{.csv.zip} archive.}
  \item{directed}{Whether to build a directed graph. \code{NULL} uses the
  archive's \code{directed} graph property, defaulting to undirected.}
}
\description{
Downloads the CSV archive of a network from
\href{https://networks.skewed.de}{Netzschleuder}, the repository used by
graph-tool's \code{collection.ns}, and loads it as a graph.
}
\details{
The download is not cached. Vertex indices in the archive are zero-based; the
returned graph, node indices, and edge endpoints are one-based. Edge
attributes, such as weights, are returned but not used by the graph, which is a
multigraph of the listed edges. Cite the original data sources given in
\code{properties} and on the dataset's Netzschleuder page.
}
\value{
A list with \code{graph}; \code{nodes}, a data frame of node properties;
\code{edges}, a data frame with \code{source}, \code{target} and any edge
properties; \code{positions}, a two-column layout matrix when the archive
provides one (otherwise \code{NULL}); and \code{properties}, a named list of
graph properties.
}
\examples{
\dontrun{
celegans <- netzschleuder_network("celegansneural")
celegans$graph
}
}
