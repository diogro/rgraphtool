\name{cpp_extensions}
\alias{cpp_extensions}
\alias{kcore_decomposition}
\alias{vertex_strength}
\alias{rgraphtool_include_dir}
\alias{rgraphtool_extension_example}
\alias{compile_cpp_extension}
\title{Native C++ extension interface, k-core decomposition, and property map access}
\usage{
kcore_decomposition(graph, deg = c("total", "out", "in"), checked = FALSE)

vertex_strength(graph, property, mode = c("total", "out", "in"),
  checked = FALSE)

rgraphtool_include_dir()

rgraphtool_extension_example(name = c("kcore", "weighted_degree"))

compile_cpp_extension(file = NULL, code = NULL, env = parent.frame(),
  rebuild = FALSE, verbose = FALSE, ...)
}
\arguments{
  \item{graph}{A graph created by \code{\link{graph_create}}.}
  \item{deg}{Degree definition for \code{kcore_decomposition()}:
  \code{"total"} (default, matching graph-tool's \code{kcore_decomposition()}),
  or \code{"out"}/\code{"in"} for directed graphs.}
  \item{property}{Character scalar naming an edge property attached to
  \code{graph} via \code{\link{graph_add_edges}} or
  \code{\link{graph_set_edge_property}}.}
  \item{mode}{Edge incidence mode for \code{vertex_strength()}: \code{"total"},
  \code{"out"}, or \code{"in"} (undirected graphs support only \code{"total"}).}
  \item{checked}{Logical; if \code{TRUE}, perform per-access bounds checking on
  vertex and edge property maps in C++; if \code{FALSE} (default), validate
  property map lengths once via \code{get_unchecked()} and use zero-overhead
  direct indexing inside hot loops.}
  \item{name}{Bundled C++ extension example name: \code{"kcore"} or
  \code{"weighted_degree"}.}
  \item{file}{Path to a C++ source file (mutually exclusive with \code{code}).}
  \item{code}{Character vector of inline C++ source code (mutually exclusive
  with \code{file}).}
  \item{env}{Target R environment where exported C++ functions are assigned.}
  \item{rebuild}{Logical; force recompilation in \code{\link[Rcpp]{sourceCpp}}.}
  \item{verbose}{Logical; print compiler invocation details.}
  \item{...}{Additional arguments forwarded to \code{\link[Rcpp]{sourceCpp}}.}
}
\description{
Exposes \pkg{rgraphtool}'s public C++ extension header (\code{<rgraphtool.h>})
for writing high-performance custom network algorithms in C++ via
\code{\link[Rcpp]{sourceCpp}} or downstream R packages declaring
\code{LinkingTo: rgraphtool, Rcpp, BH}, together with built-in reference
implementations of \eqn{k}-core decomposition (\code{kcore_decomposition()})
and weighted vertex strength (\code{vertex_strength()}).
}
\details{
\strong{Public C++ API (\code{#include <rgraphtool.h>}, namespace \code{rgraphtool}):}
\itemize{
  \item \emph{Zero-copy graph access}: \code{get_graph(SEXP obj)} accepts either
  an \code{rgraphtool_graph} R object or its external pointer, returning
  \code{Graph&}. Inspection helpers include \code{num_vertices(g)},
  \code{num_edges(g)}, \code{is_directed(g)}, \code{source(e, g)},
  \code{target(e, g)}, \code{edge_index(e)}, \code{out_degree(v, g)},
  \code{in_degree(v, g)}, \code{total_degree(v, g)}, and
  \code{checked_vertex(v, g)}.
  \item \emph{Range-based graph iteration}: \code{vertices_range(g)},
  \code{edges_range(g)}, \code{out_edges_range(v, g)},
  \code{in_edges_range(v, g)}, \code{all_edges_range(v, g)},
  \code{out_neighbors_range(v, g)}, \code{in_neighbors_range(v, g)}, and
  \code{all_neighbors_range(v, g)}. For undirected graphs,
  \code{out_edges_range(v, g)} iterates all incident edges oriented with
  \code{source == v}.
  \item \emph{Property maps (checked vs. unchecked)}:
  \code{vertex_property_view(vec, g)} wraps an R \code{NumericVector} or
  \code{IntegerVector} zero-copy in a bounds-checked \code{VertexPropertyMap<T>};
  \code{edge_property_view(g, name)} wraps an attached edge property in a
  bounds-checked \code{EdgePropertyMap<double>} / \code{ConstEdgePropertyMap<double>};
  \code{make_vertex_property<T>(g, init)} and \code{make_edge_property<T>(g, init)}
  allocate owning property maps. Calling \code{.get_unchecked()} verifies buffer
  capacity once against \code{num_vertices(g)} / \code{num_edges(g)} and returns
  a raw-pointer unchecked accessor for tight loops (analogous to graph-tool's
  \code{PropertyMap.get_unchecked()}).
}

\strong{Build requirements:} Compiling C++ extensions requires a C++23-capable
toolchain (matching \pkg{rgraphtool}'s \code{SystemRequirements: C++23}) and
the \pkg{Rcpp} and \pkg{BH} header packages. Downstream R packages can add
\code{LinkingTo: rgraphtool, Rcpp, BH} and \code{CXX_STD = CXX23} in
\code{src/Makevars}.
}
\value{
\code{kcore_decomposition()} returns an integer vector of length \eqn{N}
containing each vertex's core number. \code{vertex_strength()} returns a
numeric vector of length \eqn{N}. \code{rgraphtool_include_dir()} and
\code{rgraphtool_extension_example()} return normalized file paths.
\code{compile_cpp_extension()} returns the result of
\code{\link[Rcpp]{sourceCpp}} invisibly.
}
\examples{
g <- graph_create(6L, directed = FALSE)
graph_add_edges(g, rbind(c(1, 2), c(2, 3), c(3, 1), c(3, 4), c(4, 5)),
                properties = list(w = c(1, 2, 3, 4, 5)))
kcore_decomposition(g)
vertex_strength(g, "w")
rgraphtool_include_dir()
rgraphtool_extension_example("kcore")
}
