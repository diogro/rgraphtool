# Load the cookbook football network

Load the 115-team, 613-edge network and its supplied layout from the
graph-tool cookbook's network collection.

## Usage

``` r
football_network(file = NULL)
```

## Arguments

- file:

  Optional path to the Netzschleuder football CSV zip archive. If
  `NULL`, the archive is downloaded over HTTPS.

## Details

This function requires access to the [Netzschleuder football
dataset](https://networks.skewed.de/net/football) unless a previously
downloaded archive is supplied. The download is not cached. The source
archive stores zero-based vertex indices; returned indices and graph
vertices are one-based. The conference values are the original dataset
labels and are not inferred partitions. Node positions are the supplied
layout, not coordinates computed by this package.

The network was compiled by M. Girvan and M. E. J. Newman. Cite Girvan
and Newman (2002), doi:10.1073/pnas.122653799, when using the data.

## Value

A list with `graph`, an undirected graph; `nodes`, a data frame with
`index`, `label`, and `conference`; and `positions`, a two-column
numeric matrix in vertex order.
