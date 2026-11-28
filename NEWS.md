# rgraphtool 0.5.3

* New edge covariate type `"discrete-negbinomial"` for overdispersed counts,
  in weighted, layered and planted-partition states. It is Poisson with an
  observation-level Gamma(k, k) random effect, so counts are negative
  binomial with size `k`.
  * Microcanonical form: Dirichlet-multinomial given the block-pair total.
  * Bayesian form: a Beta prior on the NB probability.
  * At `k = 1` it equals `"discrete-geometric"` exactly.
  * On random graphs without blocks, overdispersed counts fitted as
    `"discrete-poisson"` gave 8-31 spurious blocks. Negative binomial gave 1.
  * Pass `params = list(k = )`; there is no default size.
* Weighted minimizers get weight-informed extra starts (`spectral_starts =
  TRUE`, the default) in `minimize_weighted_blockmodel_dl()`,
  `minimize_nested_weighted_blockmodel_dl()` and
  `minimize_pp_blockmodel_dl()`.
  * The search proposes merges from neighbouring groups, which carry no
    information when only the weights are structured. It stopped at one or
    two blocks even when the planted partition had a lower description length.
  * The extra starts are varimax-rotated (Rohe and Zeng 2023) and k-means
    partitions of the leading eigenvectors of the weight-centred adjacency.
  * Each is refined, and the lowest description length is kept, so fits never
    get worse; R's random state is left unchanged.
  * On graphs whose blocks show only in negative-binomial weights, NMI with the
    planted partition rose from 0.16 to 0.71, and the fits now beat the planted
    partition's description length on average.
  * Runtime roughly doubles on large graphs.
* Weighted planted-partition states have new validation tests, since
  graph-tool has no weighted PP to compare against:
  * an exact identity with graph-tool's weighted SBM on the PP buckets;
  * normalization of the weight code over integer weights;
  * detailed balance of weighted PP merge-split.
* Unused helpers no longer trigger compiler warnings (clean builds).

# rgraphtool 0.5.2

* New `pl_weight` argument for `pseudo_normal_block_state()`: a weight on the
  pseudo-log-likelihood. The default of 1 is graph-tool's model.
  * The pseudo-likelihood sums every vertex's conditional, so an edge's
    evidence is counted in both endpoint conditionals. On independent pairs
    the gain from adding an edge was 2.0 times the Gaussian likelihood gain.
  * With the quantized weight prior, extra edges cost only about 3 nats. A
    fixed share of null pairs was therefore admitted (about 7.5%), whatever
    the number of samples: 56 false edges among 36 true ones with 40 genes
    and 1000 samples.
  * `pl_weight = 0.5` is the composite-likelihood correction. It removed
    those false edges in simulations and kept full recall when samples
    outnumber genes.
  * The weight scales the data term of the description length, so compare
    description lengths only between states with the same `pl_weight`.
# rgraphtool 0.5.1

## Parallel independent starts

* New `parallel_starts()` and a `cores` argument run the `nruns`
  independent starts in forked worker processes (default
  `getOption("rgraphtool.cores", 1L)`). It is available in:
  * `compare_blockmodels()`
  * `compare_model_classes()`
  * `compare_assortative_models()`
  * `compare_weight_models()`
  * `compare_layered_models()`
  * `compare_ordered_models()`
  * the bipartite minimizers
* Each start gets its own seed drawn up front, so results are identical for
  any `cores`.
* Workers run single-threaded, and graphs are re-attached to the states
  they return. On Windows the starts run sequentially.

## Network reconstruction from an empty graph

* PseudoNormal (and other dynamics) reconstruction no longer stops at zero
  edges when edges only pay off jointly.
  * graph-tool's quantized weight prior charges a one-time entry cost when
    the graph gets its first edge: a 128-nat range constant plus
    discretization terms, about 170 nats. Any single edge gains far less, so
    the greedy search from an empty graph never added one.
  * From an empty graph, the first edge is now judged without the weight
    prior. A sweep that ends worse than where it started reverts to the
    empty graph.
  * The empty graph with fitted node parameters is kept as a fallback, so
    the result is never worse than adding no edges.
  * On a 20-node ring (M = 200), the fit went from 0 edges to the 20 true
    edges, and its DL from 5029 to 4843. On two 8-gene GMRF chains, gensbm's
    `mdl_cov` went from 0 edges to 13/14 (M = 100) and 14/14 (M = 400) true
    edges.
* Weight and node-parameter line searches expand past their +-2 window when
  the optimum lies beyond it, so unscaled data converge in fewer sweeps.

## Smaller fixes

* New S3 generic `description_length(state, ...)` for every state family,
  plus `state_blocks(state, level)`, `is_nested_state()` and `is_block_state()`, so callers need
  no class switches.
* `shrinkage_covariance()` computes its shrinkage intensity over column
  blocks. Peak memory at N = 4000 fell from 1.1 GB to 0.4 GB; results are
  unchanged (to 1e-13).
* A `parallel_for()` called from inside a parallel region now runs serially.
  Before, it deadlocked on the pool's non-recursive run lock.
* Portability:
  * The thread RNG no longer needs a 128-bit integer type.
  * Its range reduction now includes Lemire's rejection step, so it is exactly
    unbiased.
  * The vendored adjacency header no longer includes the deprecated
    `<boost/iterator.hpp>`. The package now builds without compiler warnings.
* `.gt` I/O:
  * `write_gt()` warns, naming the properties, when NA in logical or
    character properties will be written as `FALSE` / `""` (graph-tool has no
    missing values).
  * `graph_set_vertex_property()` accepts factors, as `graph_add_edges()`
    does, and stores them as character.
  * The remaining round-trip differences are documented in `?read_gt`.

## Layered models with different vertex sets per layer

* New `present` argument (a vertices x layers logical matrix) for the
  independent layers model in:
  * `layered_block_state()` and `nested_layered_block_state()`
  * the layered and nested layered minimizers
  * `compare_layered_models()`
  * `pp_block_state()` and `minimize_pp_blockmodel_dl()`
* Given presence fixes each layer's vertex set. Present vertices may have no
  edges in that layer, and presence is conditioned on rather than encoded,
  so the layer-membership term is dropped. Use it when vertex sets differ by
  design, for example genes measured in only some conditions.
* With presence inferred from the edges (the default, as in graph-tool), the
  independent-layers code over-counts. Summed over all graphs on 3 vertices,
  exp(-DL) reached 6.75. With every vertex present it sums to exactly 1.
* `layered_description_length()` now defaults to `lmemb = is.null(state$present)`.
* `present` is checked when the state is built: every edge endpoint must be
  present in the edge's layer, and named columns must be the layer levels.
* `compare_layered_models()` now errors with fewer than 2 layers.

## Weighted layered comparisons no longer depend on weight units

* The independent layers model now pools the weight hyperprior over layers,
  in `layered_block_state()`, `nested_layered_block_state()` and the
  independent `pp_block_state()`.
* Before, each layer had its own hyperprior and left its own totals
  unencoded. Multiplying all weights by s shifted the
  independent-minus-edge-covariate DL difference by:
  * (C - 1) log s for real-exponential weights;
  * 5 (C - 1) log s for flat real-normal weights.
  So the conserved-vs-rewired verdict depended on the weight units.
* The difference is now exactly invariant for flat layered and PP states.
* One fixed offset of -2 log s remains for real-normal weights when the top
  level has a single block pair: always in nested hierarchies, and in flat
  fits only at B = 1. It does not depend on the number of layers.
* Single-layer states are unchanged and still match `weighted_block_state()`.

## Bug fixes from the 0.5.1 review

* Degree description lengths:
  * All SBM families (weighted, layered, planted-partition, ranked,
    bipartite) now share one log q(n, k) table (`src/logq.h`), graph-tool's
    exact q cache plus `log_q_approx`.
  * The model-specific copies capped the exact table at 2000 edge ends and
    fell back to an incorrect approximation. Degree-corrected "distributed"
    DLs were hundreds to thousands of nats off for blocks with more than
    2000 edge ends, so the weighted / layered / PP minimizers optimized the
    wrong objective.
  * They now match graph-tool 3.9 exactly. New parity tests are in
    `test-logq-parity.R`.
  * Bipartite states no longer cast an infinite `q_max` to an integer.
* Merge-split MCMC (`multiflip = TRUE`) for weighted, planted-partition and
  layered states now satisfies detailed balance. Against exact enumeration
  of all partitions, total-variation distance dropped from 0.07-0.15 to the
  single-node noise level (0.004-0.010). Three errors were fixed:
  * Proposal probabilities: the 2^min(n, 50) cap, trivial splits treated as
    redrawn rather than rejected, a 1/C(B, 2) merge-pair probability for
    independent draws, and Gibbs-split proposal probabilities that ignored
    the random launch.
  * The number of merge-split attempts per sweep depended on the current B.
  * At finite beta the moves are now an exact uniform-split / pair-merge
    Metropolis-Hastings kernel. Greedy minimization (beta = Inf, and the PP
    minimizer) keeps the Gibbs-guided splits.
  * `sample_posterior()`, `model_evidence()` and the posterior-mode
    `compare_*()` functions were affected. New tests are in
    `test-mcmc-detailed-balance.R`.
* Nested weighted MCMC can now change the number of groups at every level
  below the top.
  * Before, singleton groups were skipped and no move created or removed a
    group, so B was fixed at every level. Posterior samples and evidence
    for nested weighted states were conditional on the starting hierarchy.
  * With probability `d` per node visit, a node now either splits off into
    a new group under its group's parent, or, if it is alone in its group,
    is absorbed into a sibling group (same parent). These two moves are
    exact reverses of each other.
  * The top level is sampled too, with one implicit parent.
  * Against exact enumeration of all 358 three-level hierarchies on 5
    vertices, TV is 0.005.
* Nested weighted and nested layered minimizers keep graph-tool's
  hierarchy depth, `ceil(log2 N) + 1` levels padded with single-group
  levels, instead of trimming to the shortest stack.
  * As in graph-tool, the number of stored levels is fixed during MCMC, and
    the effective depth adapts within the padding. With a trimmed stack the
    hierarchy could only get shallower.
  * Padding levels have zero description length. Like graph-tool's
    `print_summary()` / `get_bstack()`, `nested_summary()` and the print
    methods stop at the first single-node level and report the padding
    separately ("3 levels (+4 single-group padding)").
* PseudoNormal reconstruction:
  * The greedy sweep now optimizes edges before node parameters, as
    graph-tool does (`theta_first = False`).
  * Fitting theta to the empty graph first made every single edge look too
    costly under the quantized `xdist` prior, so reconstruction could
    return 0 edges where graph-tool recovers the network.
* `shrink_pcor()`:
  * The threaded top-k scan uses a plain C++ dot-product kernel instead of
    calling BLAS on pool workers, which oversubscribed cores with a
    multithreaded OpenBLAS. Results are bit-identical across thread counts,
    and the scan is 2x faster.
  * `nthreads` is validated.
* `nested_block_state()`: user-supplied hierarchies whose labels are not in
  first-appearance order are relabelled consistently. Before, a group could
  silently get the wrong parent.
* `pp_block_state(layers = integer(0), independent = TRUE)` on an edgeless
  graph no longer segfaults.
* `assortativity(categorical = TRUE)` compares exact values; non-integer
  values used to be truncated (0.2 and 0.7 merged). NA values are an error.
* Edge weights given to `modularity()`, `assortativity()` and other topology
  functions must be finite. NA weights, e.g. from edges added in a batch
  without the property, used to return NA silently.
* `graph_add_edges()` / `graph_set_edge_property()`:
  * Adding zero edges with properties works.
  * An R-side property on an edgeless graph keeps its length after later
    adds.
  * Factors are stored as character.
  * Changing a property between numeric and non-numeric is an error (it
    used to drop the old values silently).
* `Depends: R (>= 4.3.0)`, the first R release with C++23 support.
* `pseudo_normal_block_state()` documents `cov`, `n_eff` and `mean`.
  `shrink_pcor` is in the pkgdown index.

* Version bump marking the layered planted-partition API (`layers` /
  `independent` arguments to `pp_block_state()` and
  `minimize_pp_blockmodel_dl()`, and the `pp_edge_covariates` /
  `pp_independent_layers` options of `compare_layered_models()`), which landed
  after the 0.5.0 version number was set. Downstream packages that use the
  layered PP models should require `rgraphtool (>= 0.5.1)`.
* New `shrink_pcor()`: compiled Schafer-Strimmer shrinkage partial
  correlations (matches `corpcor::pcor.shrink()`). It uses the Woodbury
  identity in observation space, so for N >> M it avoids N x N matrices
  unless the full matrix is requested. It returns selected `pairs` or the
  `top_k` strongest pairs, and accepts a fixed `lambda` and
  `center = FALSE`. It links R's BLAS/LAPACK. The `top_k` scan is
  threaded (`nthreads`), and results do not depend on the thread count.
* `shrinkage_covariance()` gains `center = TRUE`. Use `center = FALSE` for
  observations with a known zero mean, such as regression residuals: the
  second-moment matrix `x %*% t(x) / M` is shrunk instead of the centred
  covariance.

# rgraphtool 0.5.0

Feature release adding Gram-matrix GMRF network reconstruction, multi-threaded
weighted and combined layered-weighted SBMs, and weighted planted-partition
community detection.

## Fast `PseudoNormal` Reconstruction & Merge-Split `DynBlockPrior`

* Added `shrinkage_covariance()`, a zero-dependency analytical
  Ledoit--Wolf / Schäfer--Strimmer covariance shrinkage estimator toward
  `"diagonal"` or `"identity"` targets for passing directly to
  `pseudo_normal_block_state(cov = ..., n_eff = ...)`.

* Added an exact $O(\deg)$ sufficient-statistic Gram-matrix `PseudoNormal`
  backend (`cov`, `n_eff` arguments in `pseudo_normal_block_state()`,
  auto-enabled for sample matrices when $M \ge 64, N \le 4{,}000$).
* Extended the precomputed 1D Gram edge-weight optimizer
  (`optimize_edge_weight_pn_gram()`) to both continuous Laplace and quantized
  MDL priors (`xdist = TRUE, tdist = TRUE`) via incremental $O(\log K)$
  `QuantHist` histogram maintenance.
* Parallelized candidate generation (`build_candidates()`) and added periodic
  partial-residual active-set refreshes (`refresh_every` in
  `reconstruct_dynamics()`) so greedy reconstruction recovers
  zero-marginal-correlation suppressor edges.
* Added `DynBlockPrior::merge_split_sweep()` (`merge_split = TRUE` in
  `reconstruct_dynamics()`) with BFS-seeded level-0 splits, local 2-way
  refinement, same-parent merges, and upper-level sweeps to escape single-block
  traps, plus Phase 2 finite-temperature posterior block sampling (`bs_samples`
  when `nsamples > 0`).

## Multi-Threaded Weighted & Combined Layered-Weighted SBMs

* Extracted shared conjugate/microcanonical edge-covariate likelihoods and
  sufficient statistics into `src/weight_cov.h`.
* Parallelized weighted SBM coarsening and flat refinement sweeps via
  `WeightedLevel::virtual_move_const()` and `virtual_merge_const()`, and exposed
  `nthreads` on `minimize_weighted_blockmodel_dl()` and
  `minimize_nested_weighted_blockmodel_dl()`.
* Extended `layered_block_state()`, `nested_layered_block_state()`,
  `minimize_layered_blockmodel_dl()`, `minimize_nested_layered_blockmodel_dl()`,
  `sample_posterior()`, `model_evidence()`, and `compare_layered_models()` to
  jointly model discrete edge `layers` and edge covariates (`weights`, `types`,
  `params`) under both pooled (`independent = FALSE` / `"edge_covariates"`) and
  layer-stratified (`independent = TRUE` / `"independent_layers"`) models, with
  parallel greedy sweeps (`nthreads`) and `Rcpp::Nullable` C++ default bindings.

## Weighted & Layered Planted-Partition SBMs

* Extended `pp_block_state()`, `minimize_pp_blockmodel_dl()`,
  `pp_assortativity()`, `compare_assortative_models()`, and native PP MCMC
  (`mcmc_sweep()`, `mcmc_equilibrate()`, `mcmc_anneal()`, `sample_posterior()`,
  `model_evidence()`) to support edge covariates (`weights`, `types`, `params`),
  discrete edge `layers` (`independent = FALSE` or `TRUE`), and simultaneous
  layered-weighted inference across 2 (`uniform = TRUE`) or $B + 1$
  (`uniform = FALSE`) weight buckets per stratum with $O(\deg(v))$ incremental
  moves.
* Added `pp_edge_covariates`, `pp_independent_layers`,
  `pp_edge_covariates_uniform`, and `pp_independent_layers_uniform` options to
  `compare_layered_models()`, `weight_in_mean` / `weight_out_mean` and
  per-layer (`*_layer`, `err_layer`) diagnostics to `pp_assortativity()`,
  weighted model dispatch against flat and nested weighted SBMs in
  `compare_assortative_models()`, and explicit caveats documentation in
  `?pp_block_state` (pooled inter-community `out` weights, mean-ordering / sign
  symmetry particularly at $B = 2$, and extension beyond upstream `graph-tool`
  3.9).

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
