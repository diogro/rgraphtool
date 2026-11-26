# Package index

## Graphs and data

- [`graph_create()`](https://diogro.github.io/rgraphtool/reference/graph.md)
  [`graph_add_vertices()`](https://diogro.github.io/rgraphtool/reference/graph.md)
  [`graph_add_edges()`](https://diogro.github.io/rgraphtool/reference/graph.md)
  [`graph_num_vertices()`](https://diogro.github.io/rgraphtool/reference/graph.md)
  [`graph_num_edges()`](https://diogro.github.io/rgraphtool/reference/graph.md)
  [`graph_edges()`](https://diogro.github.io/rgraphtool/reference/graph.md)
  [`graph_degree()`](https://diogro.github.io/rgraphtool/reference/graph.md)
  [`graph_set_edge_property()`](https://diogro.github.io/rgraphtool/reference/graph.md)
  [`graph_edge_property()`](https://diogro.github.io/rgraphtool/reference/graph.md)
  [`graph_edge_properties()`](https://diogro.github.io/rgraphtool/reference/graph.md)
  [`graph_set_vertex_property()`](https://diogro.github.io/rgraphtool/reference/graph.md)
  [`graph_vertex_property()`](https://diogro.github.io/rgraphtool/reference/graph.md)
  [`graph_vertex_properties()`](https://diogro.github.io/rgraphtool/reference/graph.md)
  [`graph_set_graph_property()`](https://diogro.github.io/rgraphtool/reference/graph.md)
  [`graph_graph_property()`](https://diogro.github.io/rgraphtool/reference/graph.md)
  [`graph_graph_properties()`](https://diogro.github.io/rgraphtool/reference/graph.md)
  [`print(`*`<rgraphtool_graph>`*`)`](https://diogro.github.io/rgraphtool/reference/graph.md)
  [`transform_weights()`](https://diogro.github.io/rgraphtool/reference/graph.md)
  : Create, modify, and inspect a graph and its property maps
- [`read_gt()`](https://diogro.github.io/rgraphtool/reference/graph_io.md)
  [`write_gt()`](https://diogro.github.io/rgraphtool/reference/graph_io.md)
  [`graph_load()`](https://diogro.github.io/rgraphtool/reference/graph_io.md)
  [`graph_save()`](https://diogro.github.io/rgraphtool/reference/graph_io.md)
  [`as_igraph()`](https://diogro.github.io/rgraphtool/reference/graph_io.md)
  [`from_igraph()`](https://diogro.github.io/rgraphtool/reference/graph_io.md)
  : Read, write, and convert graphs between graph-tool (.gt) and igraph
- [`football_network()`](https://diogro.github.io/rgraphtool/reference/football_network.md)
  : Load the cookbook football network
- [`fungus_tree_network()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  [`bipartite_graph()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  [`is_bipartite()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  [`as_incidence_matrix()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  [`bipartite_block_state()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  [`bipartite_description_length()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  [`minimize_bipartite_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  [`nested_bipartite_block_state()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  [`nested_bipartite_description_length()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  [`minimize_nested_bipartite_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  [`bipartite_membership()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  [`bipartite_edge_counts()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  [`bipartite_block_matrix()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  [`compare_bipartite_models()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  [`print(`*`<rgraphtool_bipartite_block_state>`*`)`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  [`print(`*`<rgraphtool_nested_bipartite_block_state>`*`)`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  [`print(`*`<rgraphtool_bipartite_model_comparison>`*`)`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  : Native bipartite network representation and Stochastic Block Model
  inference
- [`netzschleuder_network()`](https://diogro.github.io/rgraphtool/reference/netzschleuder_network.md)
  : Load a network from the Netzschleuder repository

## Block models

- [`block_state()`](https://diogro.github.io/rgraphtool/reference/block_state.md)
  [`block_membership()`](https://diogro.github.io/rgraphtool/reference/block_state.md)
  [`block_sizes()`](https://diogro.github.io/rgraphtool/reference/block_state.md)
  [`block_edge_counts()`](https://diogro.github.io/rgraphtool/reference/block_state.md)
  [`block_log_likelihood()`](https://diogro.github.io/rgraphtool/reference/block_state.md)
  [`block_description_length()`](https://diogro.github.io/rgraphtool/reference/block_state.md)
  : Inspect a fixed block-model partition
- [`fungus_tree_network()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  [`bipartite_graph()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  [`is_bipartite()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  [`as_incidence_matrix()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  [`bipartite_block_state()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  [`bipartite_description_length()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  [`minimize_bipartite_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  [`nested_bipartite_block_state()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  [`nested_bipartite_description_length()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  [`minimize_nested_bipartite_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  [`bipartite_membership()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  [`bipartite_edge_counts()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  [`bipartite_block_matrix()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  [`compare_bipartite_models()`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  [`print(`*`<rgraphtool_bipartite_block_state>`*`)`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  [`print(`*`<rgraphtool_nested_bipartite_block_state>`*`)`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  [`print(`*`<rgraphtool_bipartite_model_comparison>`*`)`](https://diogro.github.io/rgraphtool/reference/bipartite_block_state.md)
  : Native bipartite network representation and Stochastic Block Model
  inference
- [`weighted_block_state()`](https://diogro.github.io/rgraphtool/reference/weighted_block_state.md)
  [`weighted_description_length()`](https://diogro.github.io/rgraphtool/reference/weighted_block_state.md)
  [`minimize_weighted_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/weighted_block_state.md)
  [`nested_weighted_block_state()`](https://diogro.github.io/rgraphtool/reference/weighted_block_state.md)
  [`nested_weighted_level_entropy()`](https://diogro.github.io/rgraphtool/reference/weighted_block_state.md)
  [`nested_weighted_description_length()`](https://diogro.github.io/rgraphtool/reference/weighted_block_state.md)
  [`minimize_nested_weighted_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/weighted_block_state.md)
  : Flat and nested weighted stochastic block models with edge
  covariates
- [`layered_block_state()`](https://diogro.github.io/rgraphtool/reference/layered_block_state.md)
  [`layered_description_length()`](https://diogro.github.io/rgraphtool/reference/layered_block_state.md)
  [`minimize_layered_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/layered_block_state.md)
  [`nested_layered_block_state()`](https://diogro.github.io/rgraphtool/reference/layered_block_state.md)
  [`nested_layered_level_entropy()`](https://diogro.github.io/rgraphtool/reference/layered_block_state.md)
  [`nested_layered_description_length()`](https://diogro.github.io/rgraphtool/reference/layered_block_state.md)
  [`minimize_nested_layered_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/layered_block_state.md)
  : Flat and nested layered stochastic block models
- [`overlap_block_state()`](https://diogro.github.io/rgraphtool/reference/overlap_block_state.md)
  [`overlap_description_length()`](https://diogro.github.io/rgraphtool/reference/overlap_block_state.md)
  [`overlap_groups()`](https://diogro.github.io/rgraphtool/reference/overlap_block_state.md)
  [`overlap_majority_groups()`](https://diogro.github.io/rgraphtool/reference/overlap_block_state.md)
  [`minimize_overlap_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/overlap_block_state.md)
  [`nested_overlap_block_state()`](https://diogro.github.io/rgraphtool/reference/overlap_block_state.md)
  [`nested_overlap_level_entropy()`](https://diogro.github.io/rgraphtool/reference/overlap_block_state.md)
  [`nested_overlap_description_length()`](https://diogro.github.io/rgraphtool/reference/overlap_block_state.md)
  [`minimize_nested_overlap_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/overlap_block_state.md)
  [`overlap_mcmc_sweep()`](https://diogro.github.io/rgraphtool/reference/overlap_block_state.md)
  : Flat and hierarchical overlapping stochastic block models
- [`pp_block_state()`](https://diogro.github.io/rgraphtool/reference/pp_block_state.md)
  [`pp_description_length()`](https://diogro.github.io/rgraphtool/reference/pp_block_state.md)
  [`pp_assortativity()`](https://diogro.github.io/rgraphtool/reference/pp_block_state.md)
  [`minimize_pp_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/pp_block_state.md)
  [`compare_assortative_models()`](https://diogro.github.io/rgraphtool/reference/pp_block_state.md)
  : Assortative planted-partition stochastic block model
- [`ranked_block_state()`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md)
  [`ranked_description_length()`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md)
  [`minimize_ranked_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md)
  [`nested_ranked_block_state()`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md)
  [`nested_ranked_level_entropy()`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md)
  [`nested_ranked_description_length()`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md)
  [`minimize_nested_ranked_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md)
  [`ranked_vertex_order()`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md)
  [`ranked_edge_directions()`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md)
  [`ranked_edge_counts()`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md)
  [`compare_ordered_models()`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md)
  : Ordered (ranked) flat and nested stochastic block models
- [`measured_block_state()`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md)
  [`mixed_measured_block_state()`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md)
  [`uncertain_block_state()`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md)
  [`latent_multigraph_block_state()`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md)
  [`uncertain_description_length()`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md)
  [`measured_error_posterior()`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md)
  [`edge_probabilities()`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md)
  [`reconstruct_network()`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md)
  [`latent_multigraph()`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md)
  : Uncertain network reconstruction and latent Poisson multigraphs
- [`pseudo_ising_block_state()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md)
  [`ising_glauber_block_state()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md)
  [`pseudo_cising_block_state()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md)
  [`cising_glauber_block_state()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md)
  [`epidemics_block_state()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md)
  [`pseudo_normal_block_state()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md)
  [`normal_glauber_block_state()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md)
  [`linear_normal_block_state()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md)
  [`pseudo_potts_block_state()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md)
  [`potts_glauber_block_state()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md)
  [`dynamics_description_length()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md)
  [`reconstruct_dynamics()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md)
  [`reconstruct_l1_path()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md)
  [`dynamics_precision_matrix()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md)
  [`dynamics_transmission_probs()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md)
  : Network reconstruction from dynamics and behavior
- [`latent_closure_block_state()`](https://diogro.github.io/rgraphtool/reference/latent_closure_block_state.md)
  [`closure_description_length()`](https://diogro.github.io/rgraphtool/reference/latent_closure_block_state.md)
  [`infer_triadic_closure()`](https://diogro.github.io/rgraphtool/reference/latent_closure_block_state.md)
  : Latent triadic closure stochastic block model

## Inference

- [`minimize_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/minimize_blockmodel_dl.md)
  [`refine_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/minimize_blockmodel_dl.md)
  [`compare_blockmodels()`](https://diogro.github.io/rgraphtool/reference/minimize_blockmodel_dl.md)
  : Fit a flat stochastic block model by minimum description length
- [`minimize_nested_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md)
  [`nested_block_state()`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md)
  [`nested_description_length()`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md)
  [`nested_level_entropy()`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md)
  [`nested_levels()`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md)
  [`nested_project()`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md)
  [`nested_level_state()`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md)
  [`nested_summary()`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md)
  : Hierarchical (nested) stochastic block models
- [`mcmc_sweep()`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md)
  [`mcmc_equilibrate()`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md)
  [`mcmc_anneal()`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md)
  [`sample_posterior()`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md)
  : MCMC posterior sampling, equilibration, and simulated annealing for
  SBMs
- [`align_partition_labels()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md)
  [`partition_overlap()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md)
  [`align_nested_partition_labels()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md)
  [`nested_partition_overlap()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md)
  [`partition_mode_state()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md)
  [`mode_marginals()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md)
  [`mode_consensus()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md)
  [`mode_uncertainty()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md)
  [`mode_entropy()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md)
  [`mode_lprob()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md)
  [`mode_cluster_state()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md)
  [`cluster_lprob()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md)
  [`mode_summary()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md)
  : Partition label alignment, vertex marginals, consensus, and
  posterior modes
- [`predict_edges()`](https://diogro.github.io/rgraphtool/reference/predict_edges.md)
  [`evaluate_edge_prediction()`](https://diogro.github.io/rgraphtool/reference/predict_edges.md)
  : Bayesian missing and spurious edge prediction and evaluation

## Expression & preprocessing

- [`project_fixed_effects()`](https://diogro.github.io/rgraphtool/reference/expression_utils.md)
  [`filter_expression_genes()`](https://diogro.github.io/rgraphtool/reference/expression_utils.md)
  [`module_eigengenes()`](https://diogro.github.io/rgraphtool/reference/expression_utils.md)
  [`shrinkage_covariance()`](https://diogro.github.io/rgraphtool/reference/expression_utils.md)
  : Gene-Expression Fixed-Effects Projection, Filtering, Shrinkage
  Covariance, and Module Eigengenes

## Model selection

- [`compare_layered_models()`](https://diogro.github.io/rgraphtool/reference/compare_layered_models.md)
  : Compare edge-covariate and independent-layers SBM variants
- [`model_evidence()`](https://diogro.github.io/rgraphtool/reference/compare_model_classes.md)
  [`compare_model_classes()`](https://diogro.github.io/rgraphtool/reference/compare_model_classes.md)
  : Posterior model evidence and SBM model class selection
- [`compare_weight_models()`](https://diogro.github.io/rgraphtool/reference/compare_weight_models.md)
  : Compare weighted SBM covariate distributions and transformations

## Dynamics

- [`simulate_ising()`](https://diogro.github.io/rgraphtool/reference/simulate_dynamics.md)
  [`simulate_epidemics()`](https://diogro.github.io/rgraphtool/reference/simulate_dynamics.md)
  [`simulate_normal_dynamics()`](https://diogro.github.io/rgraphtool/reference/simulate_dynamics.md)
  [`simulate_potts()`](https://diogro.github.io/rgraphtool/reference/simulate_dynamics.md)
  : Native C++ forward simulators for network dynamics
- [`animate_sirs()`](https://diogro.github.io/rgraphtool/reference/animate_sirs.md)
  [`animate_dynamic_layout()`](https://diogro.github.io/rgraphtool/reference/animate_sirs.md)
  : SIRS epidemic and dynamic force-directed layout animations
- [`interactive_bfs()`](https://diogro.github.io/rgraphtool/reference/interactive_bfs.md)
  : Interactive BFS shortest-path tree explorer

## Plots and layouts

- [`plot_graph()`](https://diogro.github.io/rgraphtool/reference/plot_graph.md)
  : Draw a graph with base R graphics
- [`plot_bipartite()`](https://diogro.github.io/rgraphtool/reference/plot_bipartite.md)
  : Visualize a fitted bipartite stochastic block model
- [`project_coords()`](https://diogro.github.io/rgraphtool/reference/plot_graph_map.md)
  [`great_circle_arc()`](https://diogro.github.io/rgraphtool/reference/plot_graph_map.md)
  [`world_coastlines()`](https://diogro.github.io/rgraphtool/reference/plot_graph_map.md)
  [`plot_graph_map()`](https://diogro.github.io/rgraphtool/reference/plot_graph_map.md)
  : Geographic map projections, great-circle arcs, and spatial network
  plotting
- [`plot_hierarchy()`](https://diogro.github.io/rgraphtool/reference/plot_hierarchy.md)
  [`plot_ranked_hierarchy()`](https://diogro.github.io/rgraphtool/reference/plot_hierarchy.md)
  : Draw a nested or ranked block model hierarchy
- [`plot_marginals()`](https://diogro.github.io/rgraphtool/reference/plot_marginals.md)
  : Draw a network with pie-chart vertex marginal group probabilities
- [`prop_to_size()`](https://diogro.github.io/rgraphtool/reference/plotting_integration.md)
  [`map_colors()`](https://diogro.github.io/rgraphtool/reference/plotting_integration.md)
  [`plot_colorbar()`](https://diogro.github.io/rgraphtool/reference/plotting_integration.md)
  [`save_graph_plot()`](https://diogro.github.io/rgraphtool/reference/plotting_integration.md)
  [`graph_grob()`](https://diogro.github.io/rgraphtool/reference/plotting_integration.md)
  [`graph_to_data_frames()`](https://diogro.github.io/rgraphtool/reference/plotting_integration.md)
  [`ggplot_graph()`](https://diogro.github.io/rgraphtool/reference/plotting_integration.md)
  : Plotting integration: property scaling, colorbars, vector/raster
  export, grid, and ggplot2
- [`sfdp_layout()`](https://diogro.github.io/rgraphtool/reference/sfdp_layout.md)
  [`graph_bfs()`](https://diogro.github.io/rgraphtool/reference/sfdp_layout.md)
  : Force-directed SFDP layout and breadth-first search traversal

## Centrality, topology & generators

- [`pagerank()`](https://diogro.github.io/rgraphtool/reference/centrality.md)
  [`betweenness()`](https://diogro.github.io/rgraphtool/reference/centrality.md)
  [`closeness()`](https://diogro.github.io/rgraphtool/reference/centrality.md)
  [`eigenvector_centrality()`](https://diogro.github.io/rgraphtool/reference/centrality.md)
  [`katz_centrality()`](https://diogro.github.io/rgraphtool/reference/centrality.md)
  [`hits_centrality()`](https://diogro.github.io/rgraphtool/reference/centrality.md)
  : Vertex and edge centrality measures
- [`shortest_distance()`](https://diogro.github.io/rgraphtool/reference/topology.md)
  [`shortest_path()`](https://diogro.github.io/rgraphtool/reference/topology.md)
  [`pseudo_diameter()`](https://diogro.github.io/rgraphtool/reference/topology.md)
  [`graph_diameter()`](https://diogro.github.io/rgraphtool/reference/topology.md)
  [`label_components()`](https://diogro.github.io/rgraphtool/reference/topology.md)
  [`label_biconnected_components()`](https://diogro.github.io/rgraphtool/reference/topology.md)
  [`topological_sort()`](https://diogro.github.io/rgraphtool/reference/topology.md)
  [`is_dag()`](https://diogro.github.io/rgraphtool/reference/topology.md)
  [`local_clustering()`](https://diogro.github.io/rgraphtool/reference/topology.md)
  [`global_clustering()`](https://diogro.github.io/rgraphtool/reference/topology.md)
  [`assortativity()`](https://diogro.github.io/rgraphtool/reference/topology.md)
  : Shortest paths, connected components, clustering, and assortativity
- [`random_graph_erdos_renyi()`](https://diogro.github.io/rgraphtool/reference/random_graphs.md)
  [`random_graph_barabasi_albert()`](https://diogro.github.io/rgraphtool/reference/random_graphs.md)
  [`random_graph_watts_strogatz()`](https://diogro.github.io/rgraphtool/reference/random_graphs.md)
  [`random_graph_configuration()`](https://diogro.github.io/rgraphtool/reference/random_graphs.md)
  [`random_graph_sbm()`](https://diogro.github.io/rgraphtool/reference/random_graphs.md)
  : Classic random graph generators
- [`modularity()`](https://diogro.github.io/rgraphtool/reference/modularity.md)
  [`maximize_modularity()`](https://diogro.github.io/rgraphtool/reference/modularity.md)
  [`print(`*`<rgraphtool_modularity_state>`*`)`](https://diogro.github.io/rgraphtool/reference/modularity.md)
  : Newman modularity and modularity maximization (use with extreme
  caution)

## Package

- [`rgraphtool_max_threads()`](https://diogro.github.io/rgraphtool/reference/rgraphtool_threads.md)
  [`rgraphtool_num_threads()`](https://diogro.github.io/rgraphtool/reference/rgraphtool_threads.md)
  [`rgraphtool_set_num_threads()`](https://diogro.github.io/rgraphtool/reference/rgraphtool_threads.md)
  : Configure parallel worker threads for large-network inference
- [`kcore_decomposition()`](https://diogro.github.io/rgraphtool/reference/cpp_extensions.md)
  [`vertex_strength()`](https://diogro.github.io/rgraphtool/reference/cpp_extensions.md)
  [`rgraphtool_include_dir()`](https://diogro.github.io/rgraphtool/reference/cpp_extensions.md)
  [`rgraphtool_extension_example()`](https://diogro.github.io/rgraphtool/reference/cpp_extensions.md)
  [`compile_cpp_extension()`](https://diogro.github.io/rgraphtool/reference/cpp_extensions.md)
  : Native C++ extension interface, k-core decomposition, and property
  map access
