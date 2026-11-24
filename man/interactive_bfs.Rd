\name{interactive_bfs}
\alias{interactive_bfs}
\title{Interactive BFS shortest-path tree explorer}
\usage{
interactive_bfs(graph, pos = NULL, source = 1L, max_dist = 3L,
  directed = NULL, mode = c("html", "locator", "snapshot"), output = NULL,
  width = 760L, height = 560L, highlight_color = "#CE5C00",
  base_color = "#999999", open = FALSE,
  main = "Interactive BFS Tree Explorer")
}
\arguments{
  \item{graph}{A graph created by \code{\link{graph_create}}.}
  \item{pos}{Optional \eqn{N \times 2} vertex layout matrix; computed via
  \code{\link{sfdp_layout}} if \code{NULL}.}
  \item{source}{1-based integer vertex index (or vector of vertex indices when
  \code{mode = "snapshot"}) specifying the initial BFS root(s).}
  \item{max_dist}{Maximum BFS depth to highlight (default \code{3L}, matching
  graph-tool's \code{interactive_bst.py}).}
  \item{directed}{Logical or \code{NULL}; defaults to
  \code{isTRUE(graph$directed)}.}
  \item{mode}{Interaction mode: \code{"html"} writes a standalone zero-dependency
  HTML5 Canvas web app with live hover/click BFS tree highlighting and depth
  slider; \code{"locator"} runs an interactive base-R click loop via
  \code{\link[graphics]{locator}}; \code{"snapshot"} renders static single- or
  multi-panel BFS tree figure(s).}
  \item{output}{Optional destination HTML path (\code{mode = "html"}) or PNG
  file path (\code{mode = "snapshot"}).}
  \item{width, height}{Canvas or PNG image dimensions in pixels.}
  \item{highlight_color}{Colour for BFS tree edges and reached vertices
  (default \code{"#CE5C00"}).}
  \item{base_color}{Colour for non-tree edges and unreached vertices.}
  \item{open}{Logical; if \code{TRUE} in an interactive R session with
  \code{mode = "html"}, opens the HTML file via \code{\link[utils]{browseURL}}.}
  \item{main}{Title for the viewer or plot.}
}
\description{
R-native counterpart of graph-tool's interactive BFS visualization
(\code{interactive_bst.py}), highlighting the shortest-path BFS tree up to
depth \code{max_dist} from any selected vertex.
}
\details{
\strong{Platform and display requirements:}
\itemize{
  \item \code{mode = "html"} works headless during generation and runs
  interactively offline in any modern web browser or RStudio Viewer (hovering
  updates the BFS tree live; clicking pins or unpins the current root).
  \item \code{mode = "locator"} requires an interactive R session with a
  locator-capable graphics device (Quartz, X11, Windows, or RStudioGD); it
  errors cleanly in non-interactive sessions.
  \item \code{mode = "snapshot"} works on any R graphics device or writes a
  PNG file via \code{\link[grDevices]{png}} when \code{output} is supplied.
}
}
\value{
An object of class \code{"rgraphtool_interactive_bfs"} containing \code{bfs}
(an \code{\link{graph_bfs}} result), \code{pos}, \code{mode}, \code{output},
\code{max_dist}, and \code{directed}.
}
\examples{
g <- graph_create(10L, directed = FALSE)
graph_add_edges(g, rbind(cbind(1:9, 2:10), c(1L, 6L)))
set.seed(1)
pos <- sfdp_layout(g, max_iter = 30L)
res <- interactive_bfs(g, pos = pos, source = 1L, max_dist = 3L,
                       mode = "html", output = tempfile(fileext = ".html"))
print(res)
}
