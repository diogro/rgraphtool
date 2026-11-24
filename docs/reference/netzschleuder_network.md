# Load a network from the Netzschleuder repository

Downloads the CSV archive of a network from
[Netzschleuder](https://networks.skewed.de), the repository used by
graph-tool's `collection.ns`, and loads it as a graph.

## Usage

``` r
netzschleuder_network(name, network = name, file = NULL, directed = NULL)
```

## Arguments

- name:

  Netzschleuder dataset name, for example `"celegansneural"`.

- network:

  Network within the dataset; defaults to `name`.

- file:

  Optional path to a previously downloaded `.csv.zip` archive.

- directed:

  Whether to build a directed graph. `NULL` uses the archive's
  `directed` graph property, defaulting to undirected.

## Details

The download is not cached. Vertex indices in the archive are
zero-based; the returned graph, node indices, and edge endpoints are
one-based. Edge attributes, such as weights, are returned but not used
by the graph, which is a multigraph of the listed edges. Cite the
original data sources given in `properties` and on the dataset's
Netzschleuder page.

## Value

A list with `graph`; `nodes`, a data frame of node properties; `edges`,
a data frame with `source`, `target` and any edge properties;
`positions`, a two-column layout matrix when the archive provides one
(otherwise `NULL`); and `properties`, a named list of graph properties.

## Examples

``` r
if (FALSE) { # \dontrun{
celegans <- netzschleuder_network("celegansneural")
celegans$graph
} # }
```
