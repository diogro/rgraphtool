# rgraphtool 0.4.0

Initial public release of `rgraphtool`, a standalone native R/C++23 port of
Python `graph-tool` 3.9 with zero Python runtime dependency.

## Core Graphs, Property Maps & Interoperability

* Native `boost::adj_list` directed and undirected multigraphs with self-loops,
  parallel edges, and reference-semantics graph, vertex, and edge property maps
  (`graph_create()`, `graph_set_graph_property()`, `graph_set_vertex_property()`,
  `graph_set_edge_property()`, and `print.rgraphtool_graph()` with compact ASCII
  adjacency previews for small graphs).
* Native C++ reader and writer for `graph-tool`'s `.gt` binary format
  (`read_gt()`, `write_gt()`), supporting version 1, little- and big-endian
  encoding, all 15 property map value types, and automatic gzip/bzip2/xz/zstd
  compression.
* Lossless bidirectional `igraph` conversion (`as_igraph()`, `from_igraph()`)
  and format-detecting `graph_load()` / `graph_save()` (`.gt`, `.graphml`,
  `.gml`, `.edgelist`).
* Bundled data loaders for Netzschleuder (`netzschleuder_network()`), the
  football benchmark (`football_network()`), and the Vacher et al. (2008)
  ecological fungus--tree dataset (`fungus_tree_network()`).

## Bayesian Stochastic Block Model (SBM) Inference

* **Flat & nested (hierarchical) SBMs**: exact microcanonical description
  lengths, multilevel agglomerative minimization, and refinement
  (`block_state()`, `minimize_blockmodel_dl()`, `nested_block_state()`,
  `minimize_nested_blockmodel_dl()`, `compare_blockmodels()`).
* **Posterior MCMC & partition modes**: single-node and merge-split MCMC,
  equilibration, annealing, Hungarian label alignment, vertex marginals,
  consensus estimates, and mode clustering (`mcmc_sweep()`,
  `mcmc_equilibrate()`, `mcmc_anneal()`, `sample_posterior()`,
  `partition_mode_state()`, `mode_cluster_state()`).
* **Overlapping SBMs & model evidence**: flat and nested overlapping SBMs and
  callback-free Bayesian marginal evidence (`overlap_block_state()`,
  `nested_overlap_block_state()`, `model_evidence()`, `compare_model_classes()`).
* **Weighted & layered SBMs**: all five Peixoto (2018) edge-covariate models
  (`real-exponential`, `real-normal`, `discrete-geometric`, `discrete-poisson`,
  `discrete-binomial`) with exact change-of-variables log-Jacobians
  (`transform_weights()`, `compare_weight_models()`), plus edge-covariate and
  independent-layers SBMs (`layered_block_state()`, `compare_layered_models()`).
* **Assortative, ordered (ranked) & bipartite SBMs**: planted-partition SBMs
  (`pp_block_state()`), flat/nested ordered SBMs (`ranked_block_state()`,
  `plot_ranked_hierarchy()`), and native microcanonical flat, nested, and
  weighted bipartite SBMs (`bipartite_graph()`, `minimize_bipartite_blockmodel_dl()`,
  `minimize_nested_bipartite_blockmodel_dl()`, `plot_bipartite()`).

## Network Reconstruction

* **Uncertain networks & triadic closure**: measured, mixed-measured,
  uncertain, latent Poisson multigraph, and multi-generation latent triadic
  closure models (`measured_block_state()`, `mixed_measured_block_state()`,
  `uncertain_block_state()`, `latent_multigraph()`, `reconstruct_network()`,
  `predict_edges()`, `latent_closure_block_state()`, `infer_triadic_closure()`).
* **Reconstruction from dynamics & behavior**: all 10 dynamical models
  (equilibrium/Glauber Ising & continuous Ising, compartmental epidemics,
  Gaussian/linear-normal, and $q$-state Potts) with L1 regularization paths,
  quantized MDL priors, and native forward simulators (`reconstruct_dynamics()`,
  `reconstruct_l1_path()`, `simulate_ising()`, `simulate_epidemics()`,
  `simulate_normal_dynamics()`, `simulate_potts()`).

## Centrality, Topology, Random Graphs & Modularity

* Multi-threaded C++23 centrality (`pagerank()`, `betweenness()`, `closeness()`,
  `eigenvector_centrality()`, `katz_centrality()`, `hits_centrality()`),
  shortest paths (`shortest_distance()`, `shortest_path()`, `pseudo_diameter()`,
  `graph_diameter()`), connectivity (`label_components()`,
  `label_biconnected_components()`, `topological_sort()`, `is_dag()`),
  clustering (`local_clustering()`, `global_clustering()`), and
  `assortativity()`.
* Fast random graph generators (`random_graph_erdos_renyi()`,
  `random_graph_barabasi_albert()`, `random_graph_watts_strogatz()`,
  `random_graph_configuration()`, `random_graph_sbm()`).
* `modularity()` and `maximize_modularity()` (with explicit statistical
  overfitting warnings mirroring `graph-tool`'s `ModularityState`).

## Visualization, Scaling & C++ Extensions

* Native Barnes--Hut SFDP layouts (`sfdp_layout()`), standalone zero-dependency
  HTML5 Canvas animations and interactive BFS explorer (`animate_sirs()`,
  `animate_dynamic_layout()`, `interactive_bfs()`), base-R/`grid`/`ggplot2`
  vector/raster plotting (`plot_graph()`, `graph_grob()`, `ggplot_graph()`,
  `save_graph_plot()`), and offline cartographic map projections
  (`plot_graph_map()`).
* Portable C++23 `std::thread` parallelism (`rgraphtool_set_num_threads()`) and
  a public `<rgraphtool.h>` C++ extension interface (`LinkingTo: rgraphtool,
  Rcpp, BH`, `compile_cpp_extension()`, `kcore_decomposition()`).
