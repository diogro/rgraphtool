\name{animate_sirs}
\alias{animate_sirs}
\alias{animate_dynamic_layout}
\title{SIRS epidemic and dynamic force-directed layout animations}
\usage{
animate_sirs(graph, pos = NULL, n_steps = 120L, x = 0.001, r = 0.1, s = 0.01,
  init_state = NULL, directed = FALSE, hide_recovered = TRUE,
  fear_state = FALSE, mode = c("html", "frames", "device", "none"),
  output = NULL, width = 760L, height = 560L, fps = 12, delay = 0.05,
  open = FALSE, main = "SIRS Epidemic Animation")

animate_dynamic_layout(graph, pos = NULL, n_steps = 80L, K = 0.5, C = 0.2,
  p = 2.0, step = 0.005, n_rewire = 10L,
  mode = c("html", "frames", "device", "none"), output = NULL, width = 760L,
  height = 560L, fps = 15, delay = 0.05, open = FALSE,
  main = "Dynamic Layout & Spatial Segregation")
}
\arguments{
  \item{graph}{A graph created by \code{\link{graph_create}}.}
  \item{pos}{Optional \eqn{N \times 2} initial vertex layout matrix. Computed
  with \code{\link{sfdp_layout}} if \code{NULL}.}
  \item{n_steps}{Number of simulation steps to run (produces
  \code{n_steps + 1} frames including step 0).}
  \item{x}{Spontaneous infection probability (\code{S -> I}) per vertex per
  step.}
  \item{r}{Recovery probability (\code{I -> R}) per vertex per step.}
  \item{s}{Resusceptibility probability (\code{R -> S}) per vertex per step.}
  \item{init_state}{Optional integer vector of length \eqn{N} with values in
  \code{\{0, 1, 2\}} (\code{0 = S}, \code{1 = I}, \code{2 = R}). Defaults to
  all susceptible.}
  \item{directed}{Logical; if \code{TRUE} and \code{graph$directed} is
  \code{TRUE}, infection spreads along directed edges \code{u -> v}; otherwise
  undirected neighbors are used.}
  \item{hide_recovered}{Logical; if \code{TRUE}, recovered vertices and their
  incident edges are hidden from the animation (matching graph-tool's
  \code{animation_sirs.py}).}
  \item{fear_state}{Logical; if \code{TRUE}, susceptible vertices adjacent to
  at least one infected neighbor are highlighted in amber (matching
  graph-tool's \code{animation_zombies.py}).}
  \item{K, C, p}{SFDP spring-electrical layout parameters for
  \code{animate_dynamic_layout()}.}
  \item{step}{Per-frame incremental SFDP displacement step size.}
  \item{n_rewire}{Number of distance-shortening edge rewiring attempts per
  frame (preserving edge count and avoiding isolated vertices, self-loops, and
  parallel edges).}
  \item{mode}{Output mode: \code{"html"} writes a self-contained zero-dependency
  interactive HTML5 Canvas viewer; \code{"frames"} writes numbered PNG files to
  a directory; \code{"device"} animates on the active R graphics device;
  \code{"none"} returns the trajectory without rendering.}
  \item{output}{Optional destination HTML file path (\code{mode = "html"}) or
  directory path (\code{mode = "frames"}). Defaults to a temporary path.}
  \item{width, height}{Canvas or PNG frame dimensions in pixels.}
  \item{fps}{Initial playback rate (frames per second) for \code{mode = "html"}.}
  \item{delay}{Pause in seconds between frames when \code{mode = "device"}.}
  \item{open}{Logical; if \code{TRUE} in an interactive R session with
  \code{mode = "html"}, opens the generated HTML file via
  \code{\link[utils]{browseURL}}.}
  \item{main}{Animation title.}
}
\description{
R-native equivalents of graph-tool's SIRS epidemic animation
(\code{animation_sirs.py}, \code{animation_zombies.py}) and dynamic
force-directed spatial-segregation animation (\code{animation_dancing.py}).
Simulation loops execute natively in C++ and respect
\code{\link[base]{set.seed}}.
}
\details{
\strong{Platform and display requirements:}
\itemize{
  \item \code{mode = "html"} requires no display server, GTK+, Cairo, or
  external R packages; the generated standalone HTML5 file runs offline in any
  modern web browser or RStudio Viewer with interactive play/pause, frame
  scrubbing, speed control, and live summary charts.
  \item \code{mode = "frames"} requires R's built-in \code{\link[grDevices]{png}}
  device (available on standard macOS, Windows, and Cairo/X11-enabled Linux
  builds) and produces numbered PNG frames suitable for external encoders such
  as \code{ffmpeg}.
  \item \code{mode = "device"} renders directly to the active interactive R
  graphics device (Quartz, X11, Windows, or RStudioGD).
}
}
\value{
\code{animate_sirs()} returns an object of class
\code{"rgraphtool_sirs_animation"} containing \code{pos}, \code{states},
\code{newly_infected}, \code{fear}, \code{counts}, \code{params}, \code{mode},
\code{output}, and \code{files}.

\code{animate_dynamic_layout()} returns an object of class
\code{"rgraphtool_dynamic_layout"} containing \code{positions}, \code{edges},
\code{mean_edge_length}, \code{rewired}, \code{params}, \code{mode},
\code{output}, and \code{files}.
}
\examples{
g <- graph_create(20L, directed = FALSE)
graph_add_edges(g, cbind(1:20, c(2:20, 1L)))
set.seed(1)
init <- integer(20L); init[1L] <- 1L
sirs <- animate_sirs(g, n_steps = 10L, init_state = init, mode = "none")
print(sirs)
dyn <- animate_dynamic_layout(g, n_steps = 10L, n_rewire = 5L, mode = "none")
print(dyn)
}
