\name{random_graph_erdos_renyi}
\alias{random_graph_erdos_renyi}
\alias{random_graph_barabasi_albert}
\alias{random_graph_watts_strogatz}
\alias{random_graph_configuration}
\alias{random_graph_sbm}
\title{Classic random graph generators}
\usage{
random_graph_erdos_renyi(n, p = NULL, m = NULL, directed = FALSE,
                         self_loops = FALSE)
random_graph_barabasi_albert(n, m = 2L, gamma = 1.0, directed = FALSE)
random_graph_watts_strogatz(n, k = 4L, p = 0.05, circular = TRUE)
random_graph_configuration(out_deg, in_deg = NULL)
random_graph_sbm(blocks, probs, directed = FALSE, self_loops = FALSE)
}
\arguments{
  \item{n}{Non-negative number of vertices.}
  \item{p}{Edge probability in \eqn{[0, 1]} for \eqn{G(n, p)} or shortcut
  rewiring probability for \code{random_graph_watts_strogatz()}.}
  \item{m}{Exact number of edges for \eqn{G(n, m)}, or edges added per new
  vertex in \code{random_graph_barabasi_albert()}.}
  \item{directed}{Logical; whether the generated graph is directed.}
  \item{self_loops}{Logical; whether self-loops are permitted.}
  \item{gamma}{Preferential-attachment degree exponent
  \eqn{w_v \propto k_v^\gamma + 1} (default \code{1.0} gives linear
  Barabasi--Albert / Price attachment).}
  \item{k}{Initial lattice degree (rounded down to an even integer
  \code{2 * floor(k / 2)}) in \code{random_graph_watts_strogatz()}.}
  \item{circular}{Logical; whether the initial lattice wraps into a ring.}
  \item{out_deg, in_deg}{Non-negative integer degree sequences for stub
  matching. Omit \code{in_deg} for an undirected configuration-model multigraph.}
  \item{blocks}{Positive integer group membership vector of length \eqn{N}.}
  \item{probs}{\eqn{B \times B} matrix of inter-/intra-block edge probabilities
  in \eqn{[0, 1]}.}
}
\description{
Fast native C++23 random network generators: Erdos--Renyi \eqn{G(n, p)} (using
Batagelj--Brandes \eqn{O(n + m)} geometric skipping) and \eqn{G(n, m)} (Floyd
sampling), Barabasi--Albert / Price preferential attachment, Watts--Strogatz
small-world lattices, stub-matching configuration models, and planted
Stochastic Block Models.
}
\value{
An \code{rgraphtool_graph} object. \code{random_graph_sbm()} also attaches the
planted group memberships as vertex property \code{"block"}.
}
\examples{
set.seed(1)
g_er <- random_graph_erdos_renyi(100L, p = 0.04)
g_ba <- random_graph_barabasi_albert(100L, m = 2L)
g_ws <- random_graph_watts_strogatz(100L, k = 4L, p = 0.1)
}
