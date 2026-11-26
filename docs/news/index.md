# Changelog

## rgraphtool 0.5.0

Feature release adding native gene-expression preprocessing and
Gram-matrix GMRF network reconstruction, multi-threaded weighted and
combined layered-weighted SBMs, and weighted planted-partition community
detection.

### Gene-Expression Network Utilities

- Added zero-dependency expression preprocessing and post-clustering
  helpers:
  [`project_fixed_effects()`](https://diogro.github.io/rgraphtool/reference/expression_utils.md)
  (with `print.rgraphtool_expression_projection()`),
  [`filter_expression_genes()`](https://diogro.github.io/rgraphtool/reference/expression_utils.md),
  [`shrinkage_covariance()`](https://diogro.github.io/rgraphtool/reference/expression_utils.md)
  (analytical Ledoit–Wolf shrinkage toward `"diagonal"` or `"identity"`
  targets), and
  [`module_eigengenes()`](https://diogro.github.io/rgraphtool/reference/expression_utils.md).

### Fast `PseudoNormal` Reconstruction & Merge-Split `DynBlockPrior`

- Added an exact \\(O(\deg)\\) sufficient-statistic Gram-matrix
  `PseudoNormal` backend (`cov`, `n_eff` arguments in
  [`pseudo_normal_block_state()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md),
  auto-enabled for sample matrices when \\(M \ge 64, N \le 4{,}000\\)).
- Extended the precomputed 1D Gram edge-weight optimizer
  (`optimize_edge_weight_pn_gram()`) to both continuous Laplace and
  quantized MDL priors (`xdist = TRUE, tdist = TRUE`) via incremental
  \\(O(\log K)\\) `QuantHist` histogram maintenance.
- Parallelized candidate generation (`build_candidates()`) and added
  periodic partial-residual active-set refreshes (`refresh_every` in
  [`reconstruct_dynamics()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md))
  so greedy reconstruction recovers zero-marginal-correlation suppressor
  edges.
- Added `DynBlockPrior::merge_split_sweep()` (`merge_split = TRUE` in
  [`reconstruct_dynamics()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md))
  with BFS-seeded level-0 splits, local 2-way refinement, same-parent
  merges, and upper-level sweeps to escape single-block traps, plus
  Phase 2 finite-temperature posterior block sampling (`bs_samples` when
  `nsamples > 0`).

### Multi-Threaded Weighted & Combined Layered-Weighted SBMs

- Extracted shared conjugate/microcanonical edge-covariate likelihoods
  and sufficient statistics into `src/weight_cov.h`.
- Parallelized weighted SBM coarsening and flat refinement sweeps via
  `WeightedLevel::virtual_move_const()` and `virtual_merge_const()`, and
  exposed `nthreads` on
  [`minimize_weighted_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/weighted_block_state.md)
  and
  [`minimize_nested_weighted_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/weighted_block_state.md).
- Extended
  [`layered_block_state()`](https://diogro.github.io/rgraphtool/reference/layered_block_state.md),
  [`nested_layered_block_state()`](https://diogro.github.io/rgraphtool/reference/layered_block_state.md),
  [`minimize_layered_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/layered_block_state.md),
  [`minimize_nested_layered_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/layered_block_state.md),
  [`sample_posterior()`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md),
  [`model_evidence()`](https://diogro.github.io/rgraphtool/reference/compare_model_classes.md),
  and
  [`compare_layered_models()`](https://diogro.github.io/rgraphtool/reference/compare_layered_models.md)
  to jointly model discrete edge `layers` and edge covariates
  (`weights`, `types`, `params`) under both pooled
  (`independent = FALSE` / `"edge_covariates"`) and layer-stratified
  (`independent = TRUE` / `"independent_layers"`) models, with parallel
  greedy sweeps (`nthreads`) and `Rcpp::Nullable` C++ default bindings.

### Weighted & Layered Planted-Partition SBMs

- Extended
  [`pp_block_state()`](https://diogro.github.io/rgraphtool/reference/pp_block_state.md),
  [`minimize_pp_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/pp_block_state.md),
  [`pp_assortativity()`](https://diogro.github.io/rgraphtool/reference/pp_block_state.md),
  [`compare_assortative_models()`](https://diogro.github.io/rgraphtool/reference/pp_block_state.md),
  and native PP MCMC
  ([`mcmc_sweep()`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md),
  [`mcmc_equilibrate()`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md),
  [`mcmc_anneal()`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md),
  [`sample_posterior()`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md),
  [`model_evidence()`](https://diogro.github.io/rgraphtool/reference/compare_model_classes.md))
  to support edge covariates (`weights`, `types`, `params`), discrete
  edge `layers` (`independent = FALSE` or `TRUE`), and simultaneous
  layered-weighted inference across 2 (`uniform = TRUE`) or \\(B + 1\\)
  (`uniform = FALSE`) weight buckets per stratum with \\(O(\deg(v))\\)
  incremental moves.
- Added `pp_edge_covariates`, `pp_independent_layers`,
  `pp_edge_covariates_uniform`, and `pp_independent_layers_uniform`
  options to
  [`compare_layered_models()`](https://diogro.github.io/rgraphtool/reference/compare_layered_models.md),
  `weight_in_mean` / `weight_out_mean` and per-layer (`*_layer`,
  `err_layer`) diagnostics to
  [`pp_assortativity()`](https://diogro.github.io/rgraphtool/reference/pp_block_state.md),
  weighted model dispatch against flat and nested weighted SBMs in
  [`compare_assortative_models()`](https://diogro.github.io/rgraphtool/reference/pp_block_state.md),
  and explicit caveats documentation in
  [`?pp_block_state`](https://diogro.github.io/rgraphtool/reference/pp_block_state.md)
  (pooled inter-community `out` weights, mean-ordering / sign symmetry
  particularly at \\(B = 2\\), and extension beyond upstream
  `graph-tool` 3.9).

## rgraphtool 0.4.0

Initial public release of `rgraphtool`, a standalone native R/C++23 port
of Python `graph-tool` 3.9 with zero Python runtime dependency.

### Core Graphs, Property Maps & Interoperability

- Native `boost::adj_list` directed and undirected multigraphs with
  self-loops, parallel edges, and reference-semantics graph, vertex, and
  edge property maps
  ([`graph_create()`](https://diogro.github.io/rgraphtool/reference/graph.md),
  [`graph_set_graph_property()`](https://diogro.github.io/rgraphtool/reference/graph.md),
  [`graph_set_vertex_property()`](https://diogro.github.io/rgraphtool/reference/graph.md),
  [`graph_set_edge_property()`](https://diogro.github.io/rgraphtool/reference/graph.md),
  and
  [`print.rgraphtool_graph()`](https://diogro.github.io/rgraphtool/reference/graph.md)
  with compact ASCII adjacency previews for small graphs).
- Native C++ reader and writer for `graph-tool`’s `.gt` binary format
  ([`read_gt()`](https://diogro.github.io/rgraphtool/reference/graph_io.md),
  [`write_gt()`](https://diogro.github.io/rgraphtool/reference/graph_io.md)),
  supporting version 1, little- and big-endian encoding, all 15 property
  map value types, and automatic gzip/bzip2/xz/zstd compression.
- Lossless bidirectional `igraph` conversion
  ([`as_igraph()`](https://diogro.github.io/rgraphtool/reference/graph_io.md),
  [`from_igraph()`](https://diogro.github.io/rgraphtool/reference/graph_io.md))
  and format-detecting
  [`graph_load()`](https://diogro.github.io/rgraphtool/reference/graph_io.md)
  /
  [`graph_save()`](https://diogro.github.io/rgraphtool/reference/graph_io.md)
  (`.gt`, `.graphml`, `.gml`, `.edgelist`).
- Bundled data loaders for Netzschleuder
  ([`netzschleuder_network()`](https://diogro.github.io/rgraphtool/reference/netzschleuder_network.md)),
  the football benchmark
  ([`football_network()`](https://diogro.github.io/rgraphtool/reference/football_network.md)),
  and the Vacher et al. (2008) ecological fungus–tree dataset
  ([`fungus_tree_network()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)).

### Bayesian Stochastic Block Model (SBM) Inference

- **Flat & nested (hierarchical) SBMs**: exact microcanonical
  description lengths, multilevel agglomerative minimization, and
  refinement
  ([`block_state()`](https://diogro.github.io/rgraphtool/reference/block_state.md),
  [`minimize_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/minimize_blockmodel_dl.md),
  [`nested_block_state()`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md),
  [`minimize_nested_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md),
  [`compare_blockmodels()`](https://diogro.github.io/rgraphtool/reference/minimize_blockmodel_dl.md)).
- **Posterior MCMC & partition modes**: single-node and merge-split
  MCMC, equilibration, annealing, Hungarian label alignment, vertex
  marginals, consensus estimates, and mode clustering
  ([`mcmc_sweep()`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md),
  [`mcmc_equilibrate()`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md),
  [`mcmc_anneal()`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md),
  [`sample_posterior()`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md),
  [`partition_mode_state()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md),
  [`mode_cluster_state()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md)).
- **Overlapping SBMs & model evidence**: flat and nested overlapping
  SBMs and callback-free Bayesian marginal evidence
  ([`overlap_block_state()`](https://diogro.github.io/rgraphtool/reference/overlap_block_state.md),
  [`nested_overlap_block_state()`](https://diogro.github.io/rgraphtool/reference/overlap_block_state.md),
  [`model_evidence()`](https://diogro.github.io/rgraphtool/reference/compare_model_classes.md),
  [`compare_model_classes()`](https://diogro.github.io/rgraphtool/reference/compare_model_classes.md)).
- **Weighted & layered SBMs**: all five Peixoto (2018) edge-covariate
  models (`real-exponential`, `real-normal`, `discrete-geometric`,
  `discrete-poisson`, `discrete-binomial`) with exact
  change-of-variables log-Jacobians
  ([`transform_weights()`](https://diogro.github.io/rgraphtool/reference/graph.md),
  [`compare_weight_models()`](https://diogro.github.io/rgraphtool/reference/compare_weight_models.md)),
  plus edge-covariate and independent-layers SBMs
  ([`layered_block_state()`](https://diogro.github.io/rgraphtool/reference/layered_block_state.md),
  [`compare_layered_models()`](https://diogro.github.io/rgraphtool/reference/compare_layered_models.md)).
- **Assortative, ordered (ranked) & bipartite SBMs**: planted-partition
  SBMs
  ([`pp_block_state()`](https://diogro.github.io/rgraphtool/reference/pp_block_state.md)),
  flat/nested ordered SBMs
  ([`ranked_block_state()`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md),
  [`plot_ranked_hierarchy()`](https://diogro.github.io/rgraphtool/reference/plot_hierarchy.md)),
  and native microcanonical flat, nested, and weighted bipartite SBMs
  ([`bipartite_graph()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md),
  [`minimize_bipartite_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md),
  [`minimize_nested_bipartite_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md),
  [`plot_bipartite()`](https://diogro.github.io/rgraphtool/reference/plot_bipartite.md)).

### Network Reconstruction

- **Uncertain networks & triadic closure**: measured, mixed-measured,
  uncertain, latent Poisson multigraph, and multi-generation latent
  triadic closure models
  ([`measured_block_state()`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md),
  [`mixed_measured_block_state()`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md),
  [`uncertain_block_state()`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md),
  [`latent_multigraph()`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md),
  [`reconstruct_network()`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md),
  [`predict_edges()`](https://diogro.github.io/rgraphtool/reference/predict_edges.md),
  [`latent_closure_block_state()`](https://diogro.github.io/rgraphtool/reference/latent_closure_block_state.md),
  [`infer_triadic_closure()`](https://diogro.github.io/rgraphtool/reference/latent_closure_block_state.md)).
- **Reconstruction from dynamics & behavior**: all 10 dynamical models
  (equilibrium/Glauber Ising & continuous Ising, compartmental
  epidemics, Gaussian/linear-normal, and \\(q\\)-state Potts) with L1
  regularization paths, quantized MDL priors, and native forward
  simulators
  ([`reconstruct_dynamics()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md),
  [`reconstruct_l1_path()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md),
  [`simulate_ising()`](https://diogro.github.io/rgraphtool/reference/simulate_dynamics.md),
  [`simulate_epidemics()`](https://diogro.github.io/rgraphtool/reference/simulate_dynamics.md),
  [`simulate_normal_dynamics()`](https://diogro.github.io/rgraphtool/reference/simulate_dynamics.md),
  [`simulate_potts()`](https://diogro.github.io/rgraphtool/reference/simulate_dynamics.md)).

### Centrality, Topology, Random Graphs & Modularity

- Multi-threaded C++23 centrality
  ([`pagerank()`](https://diogro.github.io/rgraphtool/reference/centrality.md),
  [`betweenness()`](https://diogro.github.io/rgraphtool/reference/centrality.md),
  [`closeness()`](https://diogro.github.io/rgraphtool/reference/centrality.md),
  [`eigenvector_centrality()`](https://diogro.github.io/rgraphtool/reference/centrality.md),
  [`katz_centrality()`](https://diogro.github.io/rgraphtool/reference/centrality.md),
  [`hits_centrality()`](https://diogro.github.io/rgraphtool/reference/centrality.md)),
  shortest paths
  ([`shortest_distance()`](https://diogro.github.io/rgraphtool/reference/topology.md),
  [`shortest_path()`](https://diogro.github.io/rgraphtool/reference/topology.md),
  [`pseudo_diameter()`](https://diogro.github.io/rgraphtool/reference/topology.md),
  [`graph_diameter()`](https://diogro.github.io/rgraphtool/reference/topology.md)),
  connectivity
  ([`label_components()`](https://diogro.github.io/rgraphtool/reference/topology.md),
  [`label_biconnected_components()`](https://diogro.github.io/rgraphtool/reference/topology.md),
  [`topological_sort()`](https://diogro.github.io/rgraphtool/reference/topology.md),
  [`is_dag()`](https://diogro.github.io/rgraphtool/reference/topology.md)),
  clustering
  ([`local_clustering()`](https://diogro.github.io/rgraphtool/reference/topology.md),
  [`global_clustering()`](https://diogro.github.io/rgraphtool/reference/topology.md)),
  and
  [`assortativity()`](https://diogro.github.io/rgraphtool/reference/topology.md).
- Fast random graph generators
  ([`random_graph_erdos_renyi()`](https://diogro.github.io/rgraphtool/reference/random_graphs.md),
  [`random_graph_barabasi_albert()`](https://diogro.github.io/rgraphtool/reference/random_graphs.md),
  [`random_graph_watts_strogatz()`](https://diogro.github.io/rgraphtool/reference/random_graphs.md),
  [`random_graph_configuration()`](https://diogro.github.io/rgraphtool/reference/random_graphs.md),
  [`random_graph_sbm()`](https://diogro.github.io/rgraphtool/reference/random_graphs.md)).
- [`modularity()`](https://diogro.github.io/rgraphtool/reference/modularity.md)
  and
  [`maximize_modularity()`](https://diogro.github.io/rgraphtool/reference/modularity.md)
  (with explicit statistical overfitting warnings mirroring
  `graph-tool`’s `ModularityState`).

### Visualization, Scaling & C++ Extensions

- Native Barnes–Hut SFDP layouts
  ([`sfdp_layout()`](https://diogro.github.io/rgraphtool/reference/sfdp_layout.md)),
  standalone zero-dependency HTML5 Canvas animations and interactive BFS
  explorer
  ([`animate_sirs()`](https://diogro.github.io/rgraphtool/reference/animate_sirs.md),
  [`animate_dynamic_layout()`](https://diogro.github.io/rgraphtool/reference/animate_sirs.md),
  [`interactive_bfs()`](https://diogro.github.io/rgraphtool/reference/interactive_bfs.md)),
  base-R/`grid`/`ggplot2` vector/raster plotting
  ([`plot_graph()`](https://diogro.github.io/rgraphtool/reference/plot_graph.md),
  [`graph_grob()`](https://diogro.github.io/rgraphtool/reference/plotting_integration.md),
  [`ggplot_graph()`](https://diogro.github.io/rgraphtool/reference/plotting_integration.md),
  [`save_graph_plot()`](https://diogro.github.io/rgraphtool/reference/plotting_integration.md)),
  and offline cartographic map projections
  ([`plot_graph_map()`](https://diogro.github.io/rgraphtool/reference/plot_graph_map.md)).
- Portable C++23 `std::thread` parallelism
  ([`rgraphtool_set_num_threads()`](https://diogro.github.io/rgraphtool/reference/rgraphtool_threads.md))
  and a public `<rgraphtool.h>` C++ extension interface
  (`LinkingTo: rgraphtool, Rcpp, BH`,
  [`compile_cpp_extension()`](https://diogro.github.io/rgraphtool/reference/cpp_extensions.md),
  [`kcore_decomposition()`](https://diogro.github.io/rgraphtool/reference/cpp_extensions.md)).
