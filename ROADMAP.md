# rgraphtool (graph-tool-R) roadmap

## Goal and scope

Build `rgraphtool`, a standalone native R/C++ port of the Python `graph-tool`
library that can run R equivalents of the
[graph-tool cookbook](https://graph-tool.skewed.de/static/docs/stable/demos/index.html),
following its order from simpler to more advanced examples. Start with block
models and the inference guide. Python is not a runtime dependency.

The target is graph-tool 3.9 from the
[official upstream](https://git.skewed.de/count0/graph-tool), commit
`b24bc428e8426d0acb36b2d57bd61633335a3a26`. The active adjacency headers come
from this release. The flat SBM entropy and multilevel search are ports of the
3.9 inference headers (`blockmodel/`, `loops/multilevel.hh`,
`support/bisection_search.hh`); reference values come from the conda-forge
3.9 build, run in the pixi environment under `tools/reference/`. Preserve
upstream license notices and attribution when adapting code or distributing
data.

R examples should express equivalent functionality, rather than copying
Python syntax. Clearly document differences in indexing, mutable state,
defaults, output formats, and supported graph types. The matplotlib and C++
extension chapters need R-native equivalents, not a Python compatibility
layer.

**User-facing API over callbacks:** graph-tool relies heavily on user-supplied
`callback` functions inside `mcmc_equilibrate()` to collect partitions,
marginals, histograms, edge probabilities, and model-evidence traces. In
`graph-tool-R`, prefer declarative function arguments and dedicated high-level
workflows (such as `sample_posterior()`, `collect = TRUE`, `history = TRUE`,
and future `predict_edges()` / `reconstruct_network()` helpers) so cookbook
analyses run natively in C++ with ordinary arguments rather than mutable
closure callbacks. Keep `callback` minimal—an optional power-user escape hatch,
never the default path for core cookbook tasks.

## Current state

- [x] Native directed and undirected multigraph construction and inspection.
- [x] Fixed-partition block sizes, edge-count matrix, and conditional
  microcanonical likelihood.
- [x] Native flat block-model description length, with incremental move costs.
- [x] Native multilevel agglomerative minimization and refinement.
- [x] Football cookbook data loader and runnable fit/plot demonstration.
- [x] Native nested block model description length and level-by-level
  minimization, with a Netzschleuder loader, the C. elegans demo, and a
  radial hierarchy plot.
- [x] Native flat and nested single-node and merge-split MCMC, equilibration,
  annealing, Hungarian label alignment, vertex marginals, consensus estimates,
  and posterior mode clustering, with Les Misérables and hierarchical sampling
  demos and pie-chart marginal plots.
- [ ] Model class selection (overlapping/planted/weighted/layered SBMs),
  reconstruction, and general visualization.

Fixed-partition scores match graph-tool 3.9 to a relative tolerance of 1e-10 for directed
and undirected multigraphs with loops, both degree priors, and groups whose
degree sums exceed the exact integer-partition table (the Szekeres
approximation is ported). Over 20 seeds on the football network, graph-tool 3.9
finds 1738.87-1747.64 nats (mean 1740.71) without degree correction and
1787.36-1794.09 (mean 1789.18) with it. The R port finds 1738.87-1746.92
(mean 1740.08) and 1787.36-1803.30 (mean 1789.89), in 0.004-0.03 s per fit. On a
3000-vertex, 15000-edge planted partition both reach the same optimum; the R port
takes about 0.8 s against 0.8-1.7 s for graph-tool 3.9 on one thread (0.4-0.7 s
with 10 OpenMP threads).

Nested scores match graph-tool 3.9 to within 2e-11 nats on 40 fitted football
and 20 fitted C. elegans hierarchies, and on a 500-vertex multigraph with
loops, directed and undirected. Over 20 football seeds graph-tool finds
1733.53-1738.32 nats (mean 1734.23) without degree correction and
1782.02-1787.55 (mean 1782.68) with it. The port finds 1733.53-1735.20 (mean
1733.94) and 1782.02-1782.37 (mean 1782.12), in 0.01-0.02 s per fit. Over 10 C.
elegans seeds (directed), graph-tool finds 8459.65-8615.34 (mean 8506.34) and
8539.59-8662.94 (mean 8605.02). The port finds 8439.07-8578.20 (mean 8484.94)
and 8515.98-8632.54 (mean 8573.18), in 0.16-0.22 s per fit (versus 0.22-0.33 s
for single-threaded graph-tool 3.9). The sequential, lazily re-evaluated merges
are a plausible reason for the port's slightly better optima.

## Milestones

### 1. First flat SBM example

- [x] Port and validate the flat microcanonical SBM description length,
  including partition, edge-count, and degree-sequence priors. Support both
  degree-corrected and non-degree-corrected models, directed and undirected
  graphs, parallel edges, and self-loops.
- [x] Implement scalable native description-length minimization and expose fitted block
  memberships, block sizes, edge counts, and description length. Document the
  search algorithm and its differences from upstream; do not promise a global
  optimum. The multilevel search is documented in `?minimize_blockmodel_dl`;
  it runs sequentially (no OpenMP pseudo-parallel moves) and uses graph-tool's
  fast neighbour proposal.
- [x] Provide a reproducible way to obtain and load the cookbook's football
  graph, including its layout or an explicitly documented R equivalent.
- [x] Provide an R example that fits the football graph, inspects memberships
  and block counts, and produces the partition and affinity-matrix plots.

**Done when:** the example runs end to end; fixed-partition scores agree with
the pinned upstream reference within stated numerical tolerances; tests cover
edge cases and compare inferred objective values without requiring identical
stochastic labels or partitions. **Complete.** Direct group-level merge
evaluation, flat hash maps (`boost::unordered_flat_map`), dense delta buffers,
and integer `lgamma`/`log` tables bring single-threaded runtimes on par with or
faster than single-threaded graph-tool 3.9.

### 2. Hierarchical partitions and refinement

- [x] Implement nested block states, hierarchical description length,
  minimization, level traversal, and hierarchy summaries
  (`minimize_nested_blockmodel_dl()`, `nested_block_state()`,
  `nested_summary()`, `nested_project()`, `nested_level_state()`).
- [x] Load the C. elegans neural network and produce the hierarchy visualization
  and per-level membership examples (`netzschleuder_network()`,
  `plot_hierarchy()`, `demo("celegans")`). The plot is an R-native radial
  hierarchy with bundled edges, not a reproduction of graph-tool's layout.
- [x] Port the flat multilevel refinement loop (`refine=True` sweeps with
  strict merges), seeded by R's RNG.
- [x] Port the nested zero-temperature refinement loop.
- [x] Port nested merge-split MCMC and annealing at finite temperature
  (`mcmc_sweep()`, `mcmc_equilibrate()`, `mcmc_anneal()`).

**Done when:** the hierarchical cookbook examples run and refinement tests
verify that zero-temperature accepted moves do not worsen the objective.
**Complete.**

### 3. Model selection

- [x] Support comparing degree-corrected and non-degree-corrected nested models,
  including description-length differences and log posterior odds
  (`compare_blockmodels(nested = TRUE)`, `log_posterior_odds`
  \(= \Sigma_{DC} - \Sigma_{NDC} = \ln\Lambda\)). `nruns` keeps the best of
  several fits.
- [x] Run the political blogs and football model-selection examples
  (`demo("model_selection")`).

Over 5 seeds on political blogs (1490 vertices, 19090 directed edges),
graph-tool 3.9 finds 63037.5-63285.8 nats (mean 63148.4) without degree
correction and 61347.5-61644.1 (mean 61469.2) with it
(`tools/reference/polblogs_nested_upstream.csv`). The port finds
62980.2-63107.1 (mean 63043.4) and 61436.5-61558.8 (mean 61484.2), in about
1.3-1.7 s per fit (versus 1.5-3.2 s for single-threaded graph-tool 3.9). Both
clearly prefer the degree-corrected model
(\(\ln\Lambda \approx -1560\) to \(-1680\); the cookbook reports \(-1774\)).
On football both prefer the non-degree-corrected model by about 48.7 nats.

**Done when:** scores agree with the reference for identical states and the R
examples report their comparisons without assuming stochastic fits are
identical. **Complete.**

### 4. Posterior sampling

- [x] Implement the flat and nested MCMC moves and equilibration controls used
  by the guide, including merge-split sampling (`mcmc_sweep()`,
  `mcmc_equilibrate()`, `mcmc_anneal()`, `sample_posterior()`).
- [x] Add partition collection, label alignment, vertex marginals, consensus
  estimates, and posterior mode characterization (`align_partition_labels()`,
  `partition_overlap()`, `align_nested_partition_labels()`,
  `nested_partition_overlap()`, `partition_mode_state()`, `mode_marginals()`,
  `mode_consensus()`, `mode_uncertainty()`, `mode_entropy()`,
  `mode_cluster_state()`, `mode_summary()`).
- [x] Run the Les Misérables and hierarchical sampling examples with the
  corresponding visualizations (`plot_marginals()`, `demo("posterior_sampling")`,
  `demo("hierarchical_sampling")`).

Exact small-graph checks on 3-vertex directed and undirected graphs verify that
proposal probabilities normalize and match post-move reverse probabilities to
machine precision (`< 1e-12`), and that flat single-node, flat merge-split, and
nested 2-level merge-split MCMC reproduce exact Boltzmann stationary
probabilities (both uniform at \(\beta = 0\) and \(\propto e^{-\Sigma}\) at
\(\beta = 1\)). `PartitionModeState` and `ModeClusterState` description lengths,
`posterior_entropy`, `posterior_cdev`, and flat/nested `partition_overlap` match
graph-tool 3.9 to within `1e-10`. On Les Misérables (77 vertices, 254 edges),
flat equilibration and 500 thinned posterior samples (5,000 sweeps) plus mode
clustering complete in ~1.5 s; nested sampling of 400 thinned hierarchies
completes in ~0.8 s.

**Done when:** exact small-graph checks validate stationary probabilities or
transition behavior, and seeded integration examples run with sampling
uncertainty and convergence limitations documented. **Complete.**

### 5. Model class selection

- [x] Implement the model classes and model-comparison operations required by
  the cookbook's model-class selection examples, including flat and nested
  overlapping SBMs (`overlap_block_state()`, `overlap_description_length()`,
  `overlap_groups()`, `overlap_majority_groups()`,
  `minimize_overlap_blockmodel_dl()`, `nested_overlap_block_state()`,
  `nested_overlap_description_length()`, `minimize_nested_overlap_blockmodel_dl()`,
  `overlap_mcmc_sweep()`).
- [x] Add posterior partition log-probabilities (`mode_lprob()`,
  `cluster_lprob()`), nested and cluster posterior entropy (`mode_entropy()`),
  and callback-free marginal model evidence (`model_evidence()`) plus the
  high-level 8-class comparison workflow (`compare_model_classes()`).
- [x] Run the Les Misérables model-class selection cookbook workflow comparing
  MAP description lengths and full posterior marginal evidences across all
  model classes (`demo("model_class_selection")`).

Exact flat and nested overlapping description lengths match graph-tool 3.9
(`OverlapBlockState` and `NestedBlockState(base_state=OverlapBlockState)`) to
within `4e-13` on complex directed and undirected graphs (with parallel edges,
self-loops, isolated vertices, DC/NDC, distributed/uniform priors) as well as
on Les Misérables; incremental half-edge moves match full recomputation to
`< 1e-13`. `mode_lprob()`, `cluster_lprob()`, and nested/cluster `mode_entropy()`
match graph-tool 3.9's `posterior_lprob()` and `posterior_entropy()` to `1e-12`.
On Les Misérables (77 vertices, 254 edges), `demo("model_class_selection")`
completes all 8 MAP model classes (3 restarts each) plus 6 posterior evidence
integrations (200 thinned samples each) in ~5.6 s, illustrating how MAP favors
flat/overlap (`685.1`–`688.7` nats) while full Bayesian marginal evidence
decisively selects the nested model (`-log P(A) = 670.6` nats). Note that for
overlapping states, posterior entropy in `model_evidence()` is estimated from
vertex majority-group trajectories (a conservative lower bound on full
half-edge posterior entropy).

**Done when:** overlapping and non-overlapping model classes agree with the
reference for identical states, synthetic recovery tests select the true
generative class, and the cookbook comparison workflow runs end to end.
**Complete.**

### 6. Edge weights and covariates

- [x] Add edge property maps (`graph_add_edges(..., properties = ...)`,
  `graph_set_edge_property()`, `graph_edge_property()`,
  `graph_edge_properties()`, automatic attachment of numeric edge columns in
  `netzschleuder_network()`) and all five discrete and continuous Peixoto (2018)
  covariate models (`real-exponential`, `real-normal`, `discrete-geometric`,
  `discrete-poisson`, `discrete-binomial`, multiple simultaneous covariates,
  microcanonical and Bayesian conjugate hyperpriors) in native C++
  (`src/weighted.cpp`, `weighted_block_state()`, `minimize_weighted_blockmodel_dl()`,
  `nested_weighted_block_state()`, `minimize_nested_weighted_blockmodel_dl()`).
- [x] Support weight transformations with exact change-of-variables log-Jacobians
  (`transform_weights()`), declarative weighted-model comparison
  (`compare_weight_models()` in both MAP and posterior marginal evidence modes),
  flat and nested weighted MCMC / posterior sampling / `model_evidence()` /
  `plot_hierarchy()`, and run the Les Misérables weighted cookbook workflow
  (`demo("weighted_sbm")`).

Exact flat and nested weighted description lengths match graph-tool 3.9
(`WeightedBlockState` and `NestedBlockState(base_state=WeightedBlockState)`) to
within `2.84e-14` across all 5 covariate types (microcanonical, empirical-Bayes,
and custom-Bayes; DC and NDC; directed and undirected; simultaneous
multi-covariate), and incremental vertex moves match full recomputation to
`< 2.84e-14`. Single-node weighted MCMC is `1.2x`–`3.5x` faster than
single-threaded graph-tool 3.9 on `lesmis` and `celegansneural` (`0.121 s` vs
`0.291 s` for 500 `discrete-geometric` sweeps on `celegansneural`, and `0.010 s`
vs `0.035 s` for 500 `real-exponential` sweeps on `lesmis`), and the full
`demo("weighted_sbm")` workflow completes in `~1.4 s`.

**Done when:** all five covariate models and nested weighted hierarchies agree
with graph-tool 3.9 to machine precision, synthetic weight-only recovery and
Jacobian-adjusted weight-model selection pass, and the cookbook example runs end
to end. **Complete.**

### 7. Layered networks

- [x] Add discrete layer assignments (edge properties, integer/character/factor
  vectors, and categorical edge-column support in `netzschleuder_network()`) and
  native C++ flat and nested layered block states (`src/layered.cpp`,
  `layered_block_state()`, `layered_description_length()`,
  `minimize_layered_blockmodel_dl()`, `nested_layered_block_state()`,
  `nested_layered_level_entropy()`, `nested_layered_description_length()`,
  `minimize_nested_layered_blockmodel_dl()`), supporting both the edge-covariate
  (`independent = FALSE`) and conditionally independent layers
  (`independent = TRUE`) models, finite-temperature MCMC / posterior sampling
  / `model_evidence()` / `plot_hierarchy()`, and declarative model comparison
  (`compare_layered_models()` in both MAP and posterior evidence modes).
- [x] Run the layered-network cookbook workflow comparing edge-covariate vs.
  independent-layers models (MAP and full posterior marginal evidence) on Les
  Misérables (`demo("layered_sbm")`).

Exact flat and nested layered description lengths match graph-tool 3.9
(`LayeredBlockState` and `NestedBlockState(base_state=LayeredBlockState)`) to
within `2.84e-14` across all 16 reference configurations (directed and
undirected, `independent = FALSE` and `TRUE`, DC and NDC, distributed and
uniform degree priors, parallel edges within and across layers, self-loops,
isolated vertices, and blocks absent from a layer) and to `9.10e-13` on
17-layer Les Misérables; 800 incremental vertex move checks match full
recomputation to `< 2.31e-14`. Single-node layered MCMC is `2.5x`–`6.4x` faster
than single-threaded graph-tool 3.9 on `lesmis` and `celegansneural`
(`0.217 s` vs `1.384 s` for 500 `independent = TRUE` sweeps on `celegansneural`,
and `0.017 s` vs `0.067 s` for 500 `independent = FALSE` sweeps on `lesmis`),
and `demo("layered_sbm")` completes in `~3.3 s`.

**Done when:** both layered model variants and nested hierarchies agree with
graph-tool 3.9 to machine precision, synthetic layer-only community recovery
and edge-covariate vs. independent-layers model selection pass, and the
cookbook example runs end to end. **Complete.**

### 8. Assortative and ordered structures

- [x] Implement the microcanonical planted-partition model (Zhang & Peixoto
  2020) in native C++ (`src/planted_partition.cpp`, `pp_block_state()`,
  `pp_description_length()`, `pp_assortativity()`,
  `minimize_pp_blockmodel_dl()`, `compare_assortative_models()`), supporting
  both non-uniform (`uniform = FALSE`) and uniform (`uniform = TRUE`)
  assortativity, DC and NDC, distributed and uniform degree priors, MCMC /
  posterior sampling / `model_evidence()`, and the Les Misérables assortative
  cookbook workflow (`demo("assortative_sbm")`).
- [x] Implement flat and nested ordered (ranked) stochastic block models
  (Peixoto 2022) in native C++ (`src/ranked.cpp`, `ranked_block_state()`,
  `ranked_description_length()`, `minimize_ranked_blockmodel_dl()`,
  `nested_ranked_block_state()`, `nested_ranked_level_entropy()`,
  `nested_ranked_description_length()`, `minimize_nested_ranked_blockmodel_dl()`,
  `ranked_vertex_order()`, `ranked_edge_directions()`, `ranked_edge_counts()`,
  `compare_ordered_models()`), vertical rank-ordered layouts
  (`plot_ranked_hierarchy()` plus `plot_hierarchy()` support), MCMC /
  posterior sampling / `model_evidence()`, and the *C. elegans* ordered
  community cookbook workflow (`demo("ordered_sbm")`).

Exact planted-partition description lengths match graph-tool 3.9 (`PPBlockState`)
to within `1.07e-14` across all 16 reference configurations (directed and
undirected, `uniform = FALSE` and `TRUE`, DC and NDC, distributed and uniform
degree priors, parallel edges, self-loops, isolated vertices) and to `1.14e-13`
on Les Misérables; 800 incremental vertex move checks match full recomputation
to `< 2.13e-14`. Exact flat and nested ranked description lengths match
graph-tool 3.9 (`RankedBlockState` and
`NestedBlockState(base_state=RankedBlockState)`) to within `3.55e-14` across all
24 flat and 24 nested reference configurations (`distance = FALSE/TRUE`,
`unique = FALSE/TRUE`, DC/NDC, tied and distinct ranks), with 1,200 incremental
vertex-move and rank-relabel checks matching full recomputation to `< 2.62e-14`,
and matching `celegansneural` flat and nested reference entropies to `0.0`–`1.82e-12`.
On `lesmis`, planted-partition MCMC (500 sweeps) is `4.0x` faster than
single-threaded graph-tool 3.9 (`0.004 s` vs `0.016 s`); on `celegansneural`,
planted-partition minimization is `22x` faster with `174` nats lower DL
(`0.019 s` vs `0.426 s`), flat ranked MCMC (500 sweeps) is `2.2x` faster
(`0.083 s` vs `0.181 s`), and nested ranked minimization finds a `78` nats
deeper minimum (`8583.9` vs `8662.2` nats) revealing an `8:1` feed-forward flow
(`54.0%` upstream vs `6.9%` downstream). `demo("assortative_sbm")` completes in
`~0.9 s` and `demo("ordered_sbm")` completes in `~6.6 s`.

**Done when:** planted-partition and flat/nested ranked models agree with
graph-tool 3.9 to machine precision, synthetic assortative/non-assortative and
ordered/unordered model-selection tests pass, ranked layouts render cleanly, and
both cookbook examples run end to end. **Complete.**

### 9. Large networks

- [x] Build parallel and scalable inference infrastructure for very large
  networks (`N > 10^4`, `E > 10^5`): portable header-only C++23 `std::thread`
  pool (`src/parallel.h`, working identically on macOS Apple clang without
  `libomp`, Linux GCC/Clang, and Windows), thread configuration API
  (`rgraphtool_max_threads()`, `rgraphtool_num_threads()`,
  `rgraphtool_set_num_threads()`, `options(rgraphtool.nthreads = ...)`),
  lock-free read-only parallel `merge_sweep()` candidate evaluation and
  parallel batch `mh_sweep()` screening, thread-safe `LogQ` with larger
  `MathTables` (`2M` entries), dynamic-capacity upper `BlockLevel`s with
  sparse `CompactHierarchy` snapshots (eliminating `O(L * N)` upper-level
  overhead), and exact `O(deg)` incremental upper-level updates for
  `NestedOverlap` (`src/overlap.cpp`, matching full recomputation to
  `< 3.9e-14`).
- [x] Implement the guide's strictly agglomerative search
  (`control = list(bisection = FALSE, r = 0.5)`), hierarchy-shape
  initialization (`simple_init = TRUE`, `simple_init_reverse`, `base_factor`,
  `B_min_base`, `top_factor`, `strict`), and warm-start support
  (`init_blocks`, `init_levels`) on `minimize_blockmodel_dl()`,
  `refine_blockmodel_dl()`, and `minimize_nested_blockmodel_dl()`.
- [x] Run and benchmark large-network examples (`demo("large_networks")`) across
  10-core Apple Silicon and multi-core x86_64 Linux (up to 64 threads),
  documenting hardware, peak memory, and practical scaling limits.

Benchmarks on modular sparse multigraphs (`<k> = 10`):
- **10-core Apple Silicon (macOS, 32 GB RAM)**:
  - `N = 10,000, E = 49,794`: flat default `3.88 s` (1 thr) -> `0.905 s` (8 thr,
    `4.3x`, `2.8x` faster than 10-thread graph-tool 3.9's `2.54 s` with 93 nats
    lower DL); flat strictly agglomerative (`r = 0.5`) `0.69 s` (1 thr) ->
    `0.249 s` (8 thr); nested default `7.62 s` (1 thr) -> `1.85 s` (8 thr,
    `11.1x` faster than 10-thread graph-tool 3.9's `20.60 s`); nested
    `simple_init` `1.34 s` (1 thr) -> `0.534 s` (8 thr).
  - `N = 50,000, E = 249,566`: flat agglom `7.72 s` (1 thr) -> `2.14 s` (8 thr,
    `3.6x`); nested `simple_init` `16.39 s` (1 thr) -> `5.54 s` (8 thr, `3.0x`).
  - `N = 100,000, E = 499,152`: flat agglom `5.48 s` (8 thr); nested
    `simple_init` `14.82 s` (8 thr). Peak memory footprint holding both the 50k
    and 100k graphs plus `LogQ(5000)` (~100 MB) simultaneously: `604 MB`.
- **Multi-core x86_64 Linux server (GCC 11.4, R 4.5.2, tested up to 64 threads)**:
  - `N = 50,000, E = 249,566`: flat agglom `14.53 s` (1 thr) -> `6.43 s` (4 thr)
    -> `3.46 s` (16 thr, `4.2x`) -> `3.11 s` (32 thr) -> `2.99 s` (64 thr,
    `4.9x`); nested `simple_init` `31.16 s` (1 thr) -> `13.78 s` (4 thr) ->
    `8.51 s` (16 thr, `3.7x`).
  - `N = 100,000, E = 499,152`: flat agglom `32.38 s` (1 thr) -> `9.78 s`
    (16 thr) -> `7.90 s` (64 thr, `4.1x`); nested `simple_init` `71.23 s`
    (1 thr) -> `26.70 s` (16 thr) -> `22.94 s` (64 thr, `3.1x`). Peak RSS
    holding both 50k and 100k graphs and 64 thread scratch states: `791 MB`.
- **Practical scaling & limits**: `4`–`16` threads capture most of the parallel
  speedup (`3x`–`4.3x`) before sequential priority-queue merge commits dominate
  (Amdahl's law), with incremental gains up to `32`–`64` threads at
  `N >= 100,000`. On memory-constrained machines, lowering
  `options(rgraphtool.q_cache_max = 1000)` trims ~96 MB from the partition-count
  cache.

**Done when:** portable multithreading, compact nested states, incremental
`NestedOverlap` updates, strictly agglomerative search, and `simple_init` pass
correctness and scaling benchmarks up to `N = 100,000, E = 500,000`, and
`demo("large_networks")` runs end to end. **Complete.**

### 10. Uncertain network reconstruction

- [x] Implement native C++ uncertain-network reconstruction engines
  (`src/dyn_sbm.h`, `src/uncertain.cpp`, `src/closure.cpp`) with exact dynamic
  flat and nested SBM priors supporting $O(1)$ / $O(L)$ edge multiplicity edits
  and $O(\deg(v))$ vertex block moves:
  - Measured networks with unknown Beta-prior or fixed homogeneous FN/FP error
    rates (`measured_block_state()`, `measured_error_posterior()`), matching
    `MeasuredBlockState` (Peixoto 2018).
  - Heterogeneous per-dyad Beta-Binomial error models with joint hyperparameter
    MCMC (`mixed_measured_block_state()`), matching `MixedMeasuredBlockState`
    (Peixoto 2018).
  - Extraneous edge probability estimates (`uncertain_block_state()`), matching
    `UncertainBlockState` (Peixoto 2018).
  - Latent Poisson multigraphs: fast $O(V + E)$ EM (`latent_multigraph()`) and
    full Bayesian multigraph SBM (`latent_multigraph_block_state()`), matching
    `latent_multigraph()` and `LatentMultigraphBlockState` (Peixoto 2020).
  - Latent triadic closures (`latent_closure_block_state()`,
    `closure_description_length()`, `infer_triadic_closure()`), supporting
    $L \ge 1$ generations of zero-inflated triadic closures over flat or nested
    SBM substrates and per-vertex $(t, m, c)$ transitivity statistics, matching
    `LatentClosureBlockState` (Peixoto 2022).
- [x] Add exact conditional posterior edge probabilities
  (`edge_probabilities()`), joint Bayesian reconstruction
  (`reconstruct_network()`), and the declarative callback-free edge-prediction
  and evaluation workflow (`predict_edges()`, `evaluate_edge_prediction()` with
  tie-aware ROC-AUC, PR-AUC, Brier score, and ROC/PR curves), plus runnable
  cookbook demos (`demo("uncertain_reconstruction")`, `demo("triadic_closure")`).

Exact description lengths, Beta error posteriors, conditional edge
log-probabilities, and $L = 2$ triadic-closure layer entropies and $(M, E)$
triad counts match pinned graph-tool 3.9 to `< 2.4e-12` across flat and nested,
directed and undirected, degree-corrected and non-degree-corrected
configurations; incremental edge, vertex, and closure moves match full
recomputation to `< 3.2e-14`, and the edge MCMC sampler reproduces the exact
enumerated Boltzmann stationary distribution on small graphs to within Monte
Carlo error (`< 0.012`). Note two deliberate mathematical fixes relative to
graph-tool 3.9: (1) flat SBM edge deltas include exact `edges_dl` changes so
incremental `dS` matches `S(after) - S(before)` identically, and (2)
`uncertain_block_state()` uses the exact logit `log(q_default) - log1p(-q_default)`
for `q_default > 0` (fixing upstream's `log1p(+q_default)` typo). Both cookbook
demos run end to end offline in `< 0.5 s` (`demo("uncertain_reconstruction")`
achieves repeated-trial ROC-AUC `0.99` with posterior error rates `p = 0.196,
q = 0.050` vs true `0.20, 0.05`, and single-snapshot ROC-AUC `0.91` vs `0.84`
raw; `demo("triadic_closure")` separates closure vs substrate edges with
ROC-AUC `1.00`).

**Done when:** all six uncertain/latent models agree with pinned graph-tool 3.9
to machine precision, incremental moves and exact small-graph Boltzmann
stationary distributions pass, edge-prediction evaluation and triadic-closure
recovery work end to end, and both cookbook demos run cleanly. **Complete.**

### 11. Reconstruction from dynamics and behavior

- [ ] Implement the dynamical models and inference operations used by the
  synthetic and empirical examples, including L1 regularization.
- [ ] Run the indirect reconstruction examples and document input and
  preprocessing requirements.

### 12. Animation and interactive visualization

- [ ] Provide R-native equivalents for interactive animation, SIRS epidemics,
  dynamic layouts, and interactive visualizations.
- [ ] Run the examples and document any platform or display requirements.

### 13. Plotting integration

- [ ] Provide R plotting and map-integration equivalents for the matplotlib and
  basemap examples, including supported vector and raster output.
- [ ] Run the examples without requiring Python.

### 14. C++ extensions

- [ ] Document and test an R/C++ extension interface covering graph iteration,
  property map access, and checked versus unchecked access.
- [ ] Provide native R equivalents of the extension examples.

## Completion requirements for every milestone

- Add runnable, documented R examples corresponding to the relevant cookbook
  sections, with data sources, attribution, and reproducible setup.
- Verify numerical behavior against a pinned official graph-tool reference,
  using hand-checkable cases and regression fixtures as appropriate. Record
  deliberate differences and stochastic tolerances.
- Test invalid inputs, graph direction, empty graphs, loops, and parallel edges
  where applicable; run the relevant package tests and `R CMD check`.
- Specify runtime and memory limits, random-seed behavior, and platform-specific
  dependencies. Do not silently substitute a different statistical objective.
- Mark a milestone complete only when its examples execute end to end; a
  partial API or scaffold is not cookbook parity.
