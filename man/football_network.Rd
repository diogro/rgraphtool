\name{football_network}
\alias{football_network}
\title{Load the cookbook football network}
\usage{
football_network(file = NULL)
}
\arguments{
  \item{file}{Optional path to the Netzschleuder football CSV zip archive.
  If \code{NULL}, the archive is downloaded over HTTPS.}
}
\description{Load the 115-team, 613-edge network and its supplied layout from
the graph-tool cookbook's network collection.}
\details{
This function requires access to the
\href{https://networks.skewed.de/net/football}{Netzschleuder football dataset}
unless a previously downloaded archive is supplied. The download is not cached.
The source archive stores zero-based vertex indices; returned indices and
graph vertices are one-based. The conference values are the original dataset
labels and are not inferred partitions. Node positions are the supplied layout,
not coordinates computed by this package.

The network was compiled by M. Girvan and M. E. J. Newman. Cite Girvan and
Newman (2002), doi:10.1073/pnas.122653799, when using the data.
}
\value{A list with \code{graph}, an undirected graph; \code{nodes}, a data
frame with \code{index}, \code{label}, and \code{conference}; and
\code{positions}, a two-column numeric matrix in vertex order.}
