# graph-tool-R roadmap

## Goal and scope

Build a native R/C++ port that can run R equivalents of the
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

- [ ] Implement the model classes and model-comparison operations required by
  the cookbook's model-class selection examples, including overlapping
  structure where used.
- [ ] Run those examples and reproduce their reported comparison workflow.

### 6. Edge weights and covariates

- [ ] Add edge property maps and the discrete and continuous covariate models
  needed by the guide.
- [ ] Support weighted-model selection, transformations, and posterior
  sampling; run the associated examples.

### 7. Layered networks

- [ ] Add layer assignments and layered block states, inference, and model
  comparison.
- [ ] Run the layered-network cookbook examples.

### 8. Assortative and ordered structures

- [ ] Implement the planted-partition model and the assortative example.
- [ ] Implement ordered block models, ranked hierarchy layouts, and the ordered
  community example.

### 9. Large networks

- [ ] Implement the guide's strictly agglomerative search and hierarchy-shape
  initialization, including the relevant performance controls.
- [ ] Run and benchmark the large-network examples; document hardware,
  resource use, and practical limits.

### 10. Uncertain network reconstruction

- [ ] Implement measured networks, heterogeneous error models, extraneous error
  estimates, latent Poisson multigraphs, and latent triadic closures.
- [ ] Add posterior edge probabilities and the edge-prediction classification
  workflow; run the direct reconstruction examples.

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
