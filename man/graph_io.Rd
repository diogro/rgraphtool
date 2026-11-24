\name{read_gt}
\alias{read_gt}
\alias{write_gt}
\alias{graph_load}
\alias{graph_save}
\alias{as_igraph}
\alias{from_igraph}
\title{Read, write, and convert graphs between graph-tool (.gt) and igraph}
\usage{
read_gt(file)
write_gt(graph, file, comment = NULL,
         compression = c("auto", "none", "gzip", "bzip2", "xz", "zstd"))
graph_load(file, fmt = c("auto", "gt", "graphml", "gml", "edgelist"), ...)
graph_save(graph, file, fmt = c("auto", "gt", "graphml", "gml", "edgelist"), ...)
as_igraph(graph, include_properties = TRUE)
from_igraph(graph, include_properties = TRUE)
}
\arguments{
  \item{file}{Path to the input or output graph file.}
  \item{graph}{An \code{rgraphtool_graph} object (for \code{write_gt()},
  \code{graph_save()}, and \code{as_igraph()}) or an \code{igraph} object
  (for \code{from_igraph()}).}
  \item{comment}{Optional single character string stored in the \code{.gt}
  header comment field. When \code{NULL}, a descriptive default header is
  generated.}
  \item{compression}{Compression applied when writing \code{.gt} files.
  \code{"auto"} infers compression from the file extension (\code{.gz},
  \code{.bz2}, \code{.xz}, \code{.zst}); \code{read_gt()} automatically
  detects compression from magic bytes.}
  \item{fmt}{File format for \code{graph_load()} and \code{graph_save()}.
  \code{"auto"} infers the format from the file extension. \code{"gt"} uses
  the native C++ binary reader/writer without external dependencies;
  \code{"graphml"}, \code{"gml"}, and \code{"edgelist"} delegate to
  \pkg{igraph}.}
  \item{include_properties}{Logical; whether to transfer graph, vertex, and
  edge property maps during \pkg{igraph} conversion.}
  \item{...}{Additional arguments forwarded to \code{write_gt()} (for
  \code{fmt = "gt"}) or to \code{igraph::read_graph()} /
  \code{igraph::write_graph()}.}
}
\description{
Native C++ reader and writer for the \code{graph-tool} binary \code{.gt}
file format (version 1, both little-endian and big-endian, with optional
gzip/bzip2/xz/zstd compression) and lossless bidirectional conversion with
\pkg{igraph} graph objects.
}
\details{
\code{read_gt()} and \code{write_gt()} implement the official
\code{graph-tool} binary \code{.gt} specification directly in C++ with no
Python dependency. They preserve directedness, parallel edges, self-loops,
isolated vertices, and all 15 \code{graph-tool} internal property map types
(\code{bool}, \code{int16_t}, \code{int32_t}, \code{int64_t}, \code{double},
\code{long double}, \code{string}, their \code{vector<T>} counterparts, and
raw pickled \code{python::object} payloads) across graph, vertex, and edge
property maps. Numeric edge properties loaded from \code{.gt} or imported via
\code{from_igraph()} are stored directly in C++ so weighted and layered SBM
inference functions and C++ extensions can access them immediately.

\code{as_igraph()} and \code{from_igraph()} require the optional \pkg{igraph}
package and preserve vertex counts, directedness, parallel edges, self-loops,
and graph, vertex, and edge attributes.
}
\value{
\code{read_gt()}, \code{graph_load()}, and \code{from_igraph()} return an
\code{rgraphtool_graph} object. \code{as_igraph()} returns an \code{igraph}
object. \code{write_gt()} and \code{graph_save()} invisibly return \code{file}.
}
\examples{
g <- graph_create(4L, directed = TRUE)
graph_add_edges(g, rbind(c(1, 2), c(2, 3), c(3, 1), c(2, 4)),
                properties = list(weight = c(1.2, 0.8, 2.5, 1.0)))
graph_set_vertex_property(g, "name", c("a", "b", "c", "d"))
graph_set_graph_property(g, "title", "example")

f <- tempfile(fileext = ".gt.gz")
write_gt(g, f)
g2 <- read_gt(f)
graph_edge_property(g2, "weight")
graph_vertex_property(g2, "name")
unlink(f)

if (requireNamespace("igraph", quietly = TRUE)) {
  ig <- as_igraph(g)
  g3 <- from_igraph(ig)
  graph_num_edges(g3)
}
}
