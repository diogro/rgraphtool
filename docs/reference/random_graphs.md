# Classic random graph generators

Fast native C++23 random network generators: Erdos–Renyi \\G(n, p)\\
(using Batagelj–Brandes \\O(n + m)\\ geometric skipping) and \\G(n, m)\\
(Floyd sampling), Barabasi–Albert / Price preferential attachment,
Watts–Strogatz small-world lattices, stub-matching configuration models,
and planted Stochastic Block Models.

## Usage

``` r
random_graph_erdos_renyi(n, p = NULL, m = NULL, directed = FALSE,
                         self_loops = FALSE)
random_graph_barabasi_albert(n, m = 2L, gamma = 1.0, directed = FALSE)
random_graph_watts_strogatz(n, k = 4L, p = 0.05, circular = TRUE)
random_graph_configuration(out_deg, in_deg = NULL)
random_graph_sbm(blocks, probs, directed = FALSE, self_loops = FALSE)
```

## Arguments

- n:

  Non-negative number of vertices.

- p:

  Edge probability in \\\[0, 1\]\\ for \\G(n, p)\\ or shortcut rewiring
  probability for `random_graph_watts_strogatz()`.

- m:

  Exact number of edges for \\G(n, m)\\, or edges added per new vertex
  in `random_graph_barabasi_albert()`.

- directed:

  Logical; whether the generated graph is directed.

- self_loops:

  Logical; whether self-loops are permitted.

- gamma:

  Preferential-attachment degree exponent \\w_v \propto k_v^\gamma + 1\\
  (default `1.0` gives linear Barabasi–Albert / Price attachment).

- k:

  Initial lattice degree (rounded down to an even integer
  `2 * floor(k / 2)`) in `random_graph_watts_strogatz()`.

- circular:

  Logical; whether the initial lattice wraps into a ring.

- out_deg, in_deg:

  Non-negative integer degree sequences for stub matching. Omit `in_deg`
  for an undirected configuration-model multigraph.

- blocks:

  Positive integer group membership vector of length \\N\\.

- probs:

  \\B \times B\\ matrix of inter-/intra-block edge probabilities in
  \\\[0, 1\]\\.

## Value

An `rgraphtool_graph` object. `random_graph_sbm()` also attaches the
planted group memberships as vertex property `"block"`.

## Examples

``` r
set.seed(1)
g_er <- random_graph_erdos_renyi(100L, p = 0.04)
g_ba <- random_graph_barabasi_albert(100L, m = 2L)
g_ws <- random_graph_watts_strogatz(100L, k = 4L, p = 0.1)
```
