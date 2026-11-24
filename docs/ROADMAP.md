# rgraphtool (graph-tool-R) roadmap

## Goal and scope

Build `rgraphtool`, a standalone native R/C++ port of the Python
`graph-tool` library that can run R equivalents of the [graph-tool
cookbook](https://graph-tool.skewed.de/static/docs/stable/demos/index.html),
following its order from simpler to more advanced examples. Start with
block models and the inference guide. Python is not a runtime
dependency.

The target is graph-tool 3.9 from the [official
upstream](https://git.skewed.de/count0/graph-tool), commit
`b24bc428e8426d0acb36b2d57bd61633335a3a26`. The active adjacency headers
come from this release. The flat SBM entropy and multilevel search are
ports of the 3.9 inference headers (`blockmodel/`,
`loops/multilevel.hh`, `support/bisection_search.hh`); reference values
come from the conda-forge 3.9 build, run in the pixi environment under
`tools/reference/`. Preserve upstream license notices and attribution
when adapting code or distributing data.

R examples should express equivalent functionality, rather than copying
Python syntax. Clearly document differences in indexing, mutable state,
defaults, output formats, and supported graph types. The matplotlib and
C++ extension chapters need R-native equivalents, not a Python
compatibility layer.

**User-facing API over callbacks:** graph-tool relies heavily on
user-supplied `callback` functions inside
[`mcmc_equilibrate()`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md)
to collect partitions, marginals, histograms, edge probabilities, and
model-evidence traces. In `graph-tool-R`, prefer declarative function
arguments and dedicated high-level workflows (such as
[`sample_posterior()`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md),
`collect = TRUE`, `history = TRUE`, and future
[`predict_edges()`](https://diogro.github.io/rgraphtool/reference/predict_edges.md)
/
[`reconstruct_network()`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md)
helpers) so cookbook analyses run natively in C++ with ordinary
arguments rather than mutable closure callbacks. Keep `callback`
minimal—an optional power-user escape hatch, never the default path for
core cookbook tasks.

## Current state

Native directed and undirected multigraph construction and inspection.

Fixed-partition block sizes, edge-count matrix, and conditional
microcanonical likelihood.

Native flat block-model description length, with incremental move costs.

Native multilevel agglomerative minimization and refinement.

Football cookbook data loader and runnable fit/plot demonstration.

Native nested block model description length and level-by-level
minimization, with a Netzschleuder loader, the C. elegans demo, and a
radial hierarchy plot.

Native flat and nested single-node and merge-split MCMC, equilibration,
annealing, Hungarian label alignment, vertex marginals, consensus
estimates, and posterior mode clustering, with Les Misérables and
hierarchical sampling demos and pie-chart marginal plots.

Overlapping SBMs, posterior partition probabilities, and Bayesian model
class selection
([`compare_model_classes()`](https://diogro.github.io/rgraphtool/reference/compare_model_classes.md),
[`model_evidence()`](https://diogro.github.io/rgraphtool/reference/compare_model_classes.md)).

Weighted edge-covariate SBMs (all five Peixoto 2018 covariate models,
weight transformations, and Jacobian-adjusted model comparison) and
layered SBMs (edge-covariate and conditionally independent layers).

Assortative planted-partition SBMs and flat/nested ordered (ranked) SBMs
with vertical rank-ordered layouts.

Portable C++23 `std::thread` parallelism, compact nested hierarchies,
strictly agglomerative search, and `simple_init` scaling up to
`N = 100,000`.

Uncertain network reconstruction (measured, mixed-measured, uncertain,
latent Poisson multigraph, and latent triadic closure) and indirect
network reconstruction from all 10 Ising, epidemic, Gaussian, and Potts
dynamical models with L1 and quantized MDL priors.

Native SFDP force-directed layouts, BFS traversal, standalone HTML5
animations and interactive views, base-R/`grid`/`ggplot2` vector/raster
plotting, and offline geographic map projections.

Public C++ extension header interface (`inst/include/rgraphtool.h`),
range-based iteration, checked/unchecked property maps, and live
[`Rcpp::sourceCpp()`](https://rdrr.io/pkg/Rcpp/man/sourceCpp.html)
extension examples.

All 14 roadmap milestones are **complete**.

Fixed-partition scores match graph-tool 3.9 to a relative tolerance of
1e-10 for directed and undirected multigraphs with loops, both degree
priors, and groups whose degree sums exceed the exact integer-partition
table (the Szekeres approximation is ported). Over 20 seeds on the
football network, graph-tool 3.9 finds 1738.87-1747.64 nats (mean
1740.71) without degree correction and 1787.36-1794.09 (mean 1789.18)
with it. The R port finds 1738.87-1746.92 (mean 1740.08) and
1787.36-1803.30 (mean 1789.89), in 0.004-0.03 s per fit. On a
3000-vertex, 15000-edge planted partition both reach the same optimum;
the R port takes about 0.8 s against 0.8-1.7 s for graph-tool 3.9 on one
thread (0.4-0.7 s with 10 OpenMP threads).

Nested scores match graph-tool 3.9 to within 2e-11 nats on 40 fitted
football and 20 fitted C. elegans hierarchies, and on a 500-vertex
multigraph with loops, directed and undirected. Over 20 football seeds
graph-tool finds 1733.53-1738.32 nats (mean 1734.23) without degree
correction and 1782.02-1787.55 (mean 1782.68) with it. The port finds
1733.53-1735.20 (mean 1733.94) and 1782.02-1782.37 (mean 1782.12), in
0.01-0.02 s per fit. Over 10 C. elegans seeds (directed), graph-tool
finds 8459.65-8615.34 (mean 8506.34) and 8539.59-8662.94 (mean 8605.02).
The port finds 8439.07-8578.20 (mean 8484.94) and 8515.98-8632.54 (mean
8573.18), in 0.16-0.22 s per fit (versus 0.22-0.33 s for single-threaded
graph-tool 3.9). The sequential, lazily re-evaluated merges are a
plausible reason for the port’s slightly better optima.

## Milestones

### 1. First flat SBM example

Port and validate the flat microcanonical SBM description length,
including partition, edge-count, and degree-sequence priors. Support
both degree-corrected and non-degree-corrected models, directed and
undirected graphs, parallel edges, and self-loops.

Implement scalable native description-length minimization and expose
fitted block memberships, block sizes, edge counts, and description
length. Document the search algorithm and its differences from upstream;
do not promise a global optimum. The multilevel search is documented in
[`?minimize_blockmodel_dl`](https://diogro.github.io/rgraphtool/reference/minimize_blockmodel_dl.md);
it runs sequentially (no OpenMP pseudo-parallel moves) and uses
graph-tool’s fast neighbour proposal.

Provide a reproducible way to obtain and load the cookbook’s football
graph, including its layout or an explicitly documented R equivalent.

Provide an R example that fits the football graph, inspects memberships
and block counts, and produces the partition and affinity-matrix plots.

**Done when:** the example runs end to end; fixed-partition scores agree
with the pinned upstream reference within stated numerical tolerances;
tests cover edge cases and compare inferred objective values without
requiring identical stochastic labels or partitions. **Complete.**
Direct group-level merge evaluation, flat hash maps
(`boost::unordered_flat_map`), dense delta buffers, and integer
`lgamma`/`log` tables bring single-threaded runtimes on par with or
faster than single-threaded graph-tool 3.9.

### 2. Hierarchical partitions and refinement

Implement nested block states, hierarchical description length,
minimization, level traversal, and hierarchy summaries
([`minimize_nested_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md),
[`nested_block_state()`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md),
[`nested_summary()`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md),
[`nested_project()`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md),
[`nested_level_state()`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md)).

Load the C. elegans neural network and produce the hierarchy
visualization and per-level membership examples
([`netzschleuder_network()`](https://diogro.github.io/rgraphtool/reference/netzschleuder_network.md),
[`plot_hierarchy()`](https://diogro.github.io/rgraphtool/reference/plot_hierarchy.md),
`demo("celegans")`). The plot is an R-native radial hierarchy with
bundled edges, not a reproduction of graph-tool’s layout.

Port the flat multilevel refinement loop (`refine=True` sweeps with
strict merges), seeded by R’s RNG.

Port the nested zero-temperature refinement loop.

Port nested merge-split MCMC and annealing at finite temperature
([`mcmc_sweep()`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md),
[`mcmc_equilibrate()`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md),
[`mcmc_anneal()`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md)).

**Done when:** the hierarchical cookbook examples run and refinement
tests verify that zero-temperature accepted moves do not worsen the
objective. **Complete.**

### 3. Model selection

Support comparing degree-corrected and non-degree-corrected nested
models, including description-length differences and log posterior odds
(`compare_blockmodels(nested = TRUE)`, `log_posterior_odds` (= *{DC} -*
{NDC} = )). `nruns` keeps the best of several fits.

Run the political blogs and football model-selection examples
(`demo("model_selection")`).

Over 5 seeds on political blogs (1490 vertices, 19090 directed edges),
graph-tool 3.9 finds 63037.5-63285.8 nats (mean 63148.4) without degree
correction and 61347.5-61644.1 (mean 61469.2) with it
(`tools/reference/polblogs_nested_upstream.csv`). The port finds
62980.2-63107.1 (mean 63043.4) and 61436.5-61558.8 (mean 61484.2), in
about 1.3-1.7 s per fit (versus 1.5-3.2 s for single-threaded graph-tool
3.9). Both clearly prefer the degree-corrected model (() to (-1680); the
cookbook reports (-1774)). On football both prefer the
non-degree-corrected model by about 48.7 nats.

**Done when:** scores agree with the reference for identical states and
the R examples report their comparisons without assuming stochastic fits
are identical. **Complete.**

### 4. Posterior sampling

Implement the flat and nested MCMC moves and equilibration controls used
by the guide, including merge-split sampling
([`mcmc_sweep()`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md),
[`mcmc_equilibrate()`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md),
[`mcmc_anneal()`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md),
[`sample_posterior()`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md)).

Add partition collection, label alignment, vertex marginals, consensus
estimates, and posterior mode characterization
([`align_partition_labels()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md),
[`partition_overlap()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md),
[`align_nested_partition_labels()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md),
[`nested_partition_overlap()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md),
[`partition_mode_state()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md),
[`mode_marginals()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md),
[`mode_consensus()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md),
[`mode_uncertainty()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md),
[`mode_entropy()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md),
[`mode_cluster_state()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md),
[`mode_summary()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md)).

Run the Les Misérables and hierarchical sampling examples with the
corresponding visualizations
([`plot_marginals()`](https://diogro.github.io/rgraphtool/reference/plot_marginals.md),
`demo("posterior_sampling")`, `demo("hierarchical_sampling")`).

Exact small-graph checks on 3-vertex directed and undirected graphs
verify that proposal probabilities normalize and match post-move reverse
probabilities to machine precision (`< 1e-12`), and that flat
single-node, flat merge-split, and nested 2-level merge-split MCMC
reproduce exact Boltzmann stationary probabilities (both uniform at (=
0) and (e^{-}) at (= 1)). `PartitionModeState` and `ModeClusterState`
description lengths, `posterior_entropy`, `posterior_cdev`, and
flat/nested `partition_overlap` match graph-tool 3.9 to within `1e-10`.
On Les Misérables (77 vertices, 254 edges), flat equilibration and 500
thinned posterior samples (5,000 sweeps) plus mode clustering complete
in ~1.5 s; nested sampling of 400 thinned hierarchies completes in ~0.8
s.

**Done when:** exact small-graph checks validate stationary
probabilities or transition behavior, and seeded integration examples
run with sampling uncertainty and convergence limitations documented.
**Complete.**

### 5. Model class selection

Implement the model classes and model-comparison operations required by
the cookbook’s model-class selection examples, including flat and nested
overlapping SBMs
([`overlap_block_state()`](https://diogro.github.io/rgraphtool/reference/overlap_block_state.md),
[`overlap_description_length()`](https://diogro.github.io/rgraphtool/reference/overlap_block_state.md),
[`overlap_groups()`](https://diogro.github.io/rgraphtool/reference/overlap_block_state.md),
[`overlap_majority_groups()`](https://diogro.github.io/rgraphtool/reference/overlap_block_state.md),
[`minimize_overlap_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/overlap_block_state.md),
[`nested_overlap_block_state()`](https://diogro.github.io/rgraphtool/reference/overlap_block_state.md),
[`nested_overlap_description_length()`](https://diogro.github.io/rgraphtool/reference/overlap_block_state.md),
[`minimize_nested_overlap_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/overlap_block_state.md),
[`overlap_mcmc_sweep()`](https://diogro.github.io/rgraphtool/reference/overlap_block_state.md)).

Add posterior partition log-probabilities
([`mode_lprob()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md),
[`cluster_lprob()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md)),
nested and cluster posterior entropy
([`mode_entropy()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md)),
and callback-free marginal model evidence
([`model_evidence()`](https://diogro.github.io/rgraphtool/reference/compare_model_classes.md))
plus the high-level 8-class comparison workflow
([`compare_model_classes()`](https://diogro.github.io/rgraphtool/reference/compare_model_classes.md)).

Run the Les Misérables model-class selection cookbook workflow comparing
MAP description lengths and full posterior marginal evidences across all
model classes (`demo("model_class_selection")`).

Exact flat and nested overlapping description lengths match graph-tool
3.9 (`OverlapBlockState` and
`NestedBlockState(base_state=OverlapBlockState)`) to within `4e-13` on
complex directed and undirected graphs (with parallel edges, self-loops,
isolated vertices, DC/NDC, distributed/uniform priors) as well as on Les
Misérables; incremental half-edge moves match full recomputation to
`< 1e-13`.
[`mode_lprob()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md),
[`cluster_lprob()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md),
and nested/cluster
[`mode_entropy()`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md)
match graph-tool 3.9’s `posterior_lprob()` and `posterior_entropy()` to
`1e-12`. On Les Misérables (77 vertices, 254 edges),
`demo("model_class_selection")` completes all 8 MAP model classes (3
restarts each) plus 6 posterior evidence integrations (200 thinned
samples each) in ~5.6 s, illustrating how MAP favors flat/overlap
(`685.1`–`688.7` nats) while full Bayesian marginal evidence decisively
selects the nested model (`-log P(A) = 670.6` nats). Note that for
overlapping states, posterior entropy in
[`model_evidence()`](https://diogro.github.io/rgraphtool/reference/compare_model_classes.md)
is estimated from vertex majority-group trajectories (a conservative
lower bound on full half-edge posterior entropy).

**Done when:** overlapping and non-overlapping model classes agree with
the reference for identical states, synthetic recovery tests select the
true generative class, and the cookbook comparison workflow runs end to
end. **Complete.**

### 6. Edge weights and covariates

Add edge property maps (`graph_add_edges(..., properties = ...)`,
[`graph_set_edge_property()`](https://diogro.github.io/rgraphtool/reference/graph.md),
[`graph_edge_property()`](https://diogro.github.io/rgraphtool/reference/graph.md),
[`graph_edge_properties()`](https://diogro.github.io/rgraphtool/reference/graph.md),
automatic attachment of numeric edge columns in
[`netzschleuder_network()`](https://diogro.github.io/rgraphtool/reference/netzschleuder_network.md))
and all five discrete and continuous Peixoto (2018) covariate models
(`real-exponential`, `real-normal`, `discrete-geometric`,
`discrete-poisson`, `discrete-binomial`, multiple simultaneous
covariates, microcanonical and Bayesian conjugate hyperpriors) in native
C++ (`src/weighted.cpp`,
[`weighted_block_state()`](https://diogro.github.io/rgraphtool/reference/weighted_block_state.md),
[`minimize_weighted_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/weighted_block_state.md),
[`nested_weighted_block_state()`](https://diogro.github.io/rgraphtool/reference/weighted_block_state.md),
[`minimize_nested_weighted_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/weighted_block_state.md)).

Support weight transformations with exact change-of-variables
log-Jacobians
([`transform_weights()`](https://diogro.github.io/rgraphtool/reference/graph.md)),
declarative weighted-model comparison
([`compare_weight_models()`](https://diogro.github.io/rgraphtool/reference/compare_weight_models.md)
in both MAP and posterior marginal evidence modes), flat and nested
weighted MCMC / posterior sampling /
[`model_evidence()`](https://diogro.github.io/rgraphtool/reference/compare_model_classes.md)
/
[`plot_hierarchy()`](https://diogro.github.io/rgraphtool/reference/plot_hierarchy.md),
and run the Les Misérables weighted cookbook workflow
(`demo("weighted_sbm")`).

Exact flat and nested weighted description lengths match graph-tool 3.9
(`WeightedBlockState` and
`NestedBlockState(base_state=WeightedBlockState)`) to within `2.84e-14`
across all 5 covariate types (microcanonical, empirical-Bayes, and
custom-Bayes; DC and NDC; directed and undirected; simultaneous
multi-covariate), and incremental vertex moves match full recomputation
to `< 2.84e-14`. Single-node weighted MCMC is `1.2x`–`3.5x` faster than
single-threaded graph-tool 3.9 on `lesmis` and `celegansneural`
(`0.121 s` vs `0.291 s` for 500 `discrete-geometric` sweeps on
`celegansneural`, and `0.010 s` vs `0.035 s` for 500 `real-exponential`
sweeps on `lesmis`), and the full `demo("weighted_sbm")` workflow
completes in `~1.4 s`.

**Done when:** all five covariate models and nested weighted hierarchies
agree with graph-tool 3.9 to machine precision, synthetic weight-only
recovery and Jacobian-adjusted weight-model selection pass, and the
cookbook example runs end to end. **Complete.**

### 7. Layered networks

Add discrete layer assignments (edge properties,
integer/character/factor vectors, and categorical edge-column support in
[`netzschleuder_network()`](https://diogro.github.io/rgraphtool/reference/netzschleuder_network.md))
and native C++ flat and nested layered block states (`src/layered.cpp`,
[`layered_block_state()`](https://diogro.github.io/rgraphtool/reference/layered_block_state.md),
[`layered_description_length()`](https://diogro.github.io/rgraphtool/reference/layered_block_state.md),
[`minimize_layered_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/layered_block_state.md),
[`nested_layered_block_state()`](https://diogro.github.io/rgraphtool/reference/layered_block_state.md),
[`nested_layered_level_entropy()`](https://diogro.github.io/rgraphtool/reference/layered_block_state.md),
[`nested_layered_description_length()`](https://diogro.github.io/rgraphtool/reference/layered_block_state.md),
[`minimize_nested_layered_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/layered_block_state.md)),
supporting both the edge-covariate (`independent = FALSE`) and
conditionally independent layers (`independent = TRUE`) models,
finite-temperature MCMC / posterior sampling /
[`model_evidence()`](https://diogro.github.io/rgraphtool/reference/compare_model_classes.md)
/
[`plot_hierarchy()`](https://diogro.github.io/rgraphtool/reference/plot_hierarchy.md),
and declarative model comparison
([`compare_layered_models()`](https://diogro.github.io/rgraphtool/reference/compare_layered_models.md)
in both MAP and posterior evidence modes).

Run the layered-network cookbook workflow comparing edge-covariate vs.
independent-layers models (MAP and full posterior marginal evidence) on
Les Misérables (`demo("layered_sbm")`).

Exact flat and nested layered description lengths match graph-tool 3.9
(`LayeredBlockState` and
`NestedBlockState(base_state=LayeredBlockState)`) to within `2.84e-14`
across all 16 reference configurations (directed and undirected,
`independent = FALSE` and `TRUE`, DC and NDC, distributed and uniform
degree priors, parallel edges within and across layers, self-loops,
isolated vertices, and blocks absent from a layer) and to `9.10e-13` on
17-layer Les Misérables; 800 incremental vertex move checks match full
recomputation to `< 2.31e-14`. Single-node layered MCMC is `2.5x`–`6.4x`
faster than single-threaded graph-tool 3.9 on `lesmis` and
`celegansneural` (`0.217 s` vs `1.384 s` for 500 `independent = TRUE`
sweeps on `celegansneural`, and `0.017 s` vs `0.067 s` for 500
`independent = FALSE` sweeps on `lesmis`), and `demo("layered_sbm")`
completes in `~3.3 s`.

**Done when:** both layered model variants and nested hierarchies agree
with graph-tool 3.9 to machine precision, synthetic layer-only community
recovery and edge-covariate vs. independent-layers model selection pass,
and the cookbook example runs end to end. **Complete.**

### 8. Assortative and ordered structures

Implement the microcanonical planted-partition model (Zhang & Peixoto

2020. in native C++ (`src/planted_partition.cpp`,
      [`pp_block_state()`](https://diogro.github.io/rgraphtool/reference/pp_block_state.md),
      [`pp_description_length()`](https://diogro.github.io/rgraphtool/reference/pp_block_state.md),
      [`pp_assortativity()`](https://diogro.github.io/rgraphtool/reference/pp_block_state.md),
      [`minimize_pp_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/pp_block_state.md),
      [`compare_assortative_models()`](https://diogro.github.io/rgraphtool/reference/pp_block_state.md)),
      supporting both non-uniform (`uniform = FALSE`) and uniform
      (`uniform = TRUE`) assortativity, DC and NDC, distributed and
      uniform degree priors, MCMC / posterior sampling /
      [`model_evidence()`](https://diogro.github.io/rgraphtool/reference/compare_model_classes.md),
      and the Les Misérables assortative cookbook workflow
      (`demo("assortative_sbm")`).

Implement flat and nested ordered (ranked) stochastic block models
(Peixoto 2022) in native C++ (`src/ranked.cpp`,
[`ranked_block_state()`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md),
[`ranked_description_length()`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md),
[`minimize_ranked_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md),
[`nested_ranked_block_state()`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md),
[`nested_ranked_level_entropy()`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md),
[`nested_ranked_description_length()`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md),
[`minimize_nested_ranked_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md),
[`ranked_vertex_order()`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md),
[`ranked_edge_directions()`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md),
[`ranked_edge_counts()`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md),
[`compare_ordered_models()`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md)),
vertical rank-ordered layouts
([`plot_ranked_hierarchy()`](https://diogro.github.io/rgraphtool/reference/plot_hierarchy.md)
plus
[`plot_hierarchy()`](https://diogro.github.io/rgraphtool/reference/plot_hierarchy.md)
support), MCMC / posterior sampling /
[`model_evidence()`](https://diogro.github.io/rgraphtool/reference/compare_model_classes.md),
and the *C. elegans* ordered community cookbook workflow
(`demo("ordered_sbm")`).

Exact planted-partition description lengths match graph-tool 3.9
(`PPBlockState`) to within `1.07e-14` across all 16 reference
configurations (directed and undirected, `uniform = FALSE` and `TRUE`,
DC and NDC, distributed and uniform degree priors, parallel edges,
self-loops, isolated vertices) and to `1.14e-13` on Les Misérables; 800
incremental vertex move checks match full recomputation to `< 2.13e-14`.
Exact flat and nested ranked description lengths match graph-tool 3.9
(`RankedBlockState` and `NestedBlockState(base_state=RankedBlockState)`)
to within `3.55e-14` across all 24 flat and 24 nested reference
configurations (`distance = FALSE/TRUE`, `unique = FALSE/TRUE`, DC/NDC,
tied and distinct ranks), with 1,200 incremental vertex-move and
rank-relabel checks matching full recomputation to `< 2.62e-14`, and
matching `celegansneural` flat and nested reference entropies to
`0.0`–`1.82e-12`. On `lesmis`, planted-partition MCMC (500 sweeps) is
`4.0x` faster than single-threaded graph-tool 3.9 (`0.004 s` vs
`0.016 s`); on `celegansneural`, planted-partition minimization is `22x`
faster with `174` nats lower DL (`0.019 s` vs `0.426 s`), flat ranked
MCMC (500 sweeps) is `2.2x` faster (`0.083 s` vs `0.181 s`), and nested
ranked minimization finds a `78` nats deeper minimum (`8583.9` vs
`8662.2` nats) revealing an `8:1` feed-forward flow (`54.0%` upstream vs
`6.9%` downstream). `demo("assortative_sbm")` completes in `~0.9 s` and
`demo("ordered_sbm")` completes in `~6.6 s`.

**Done when:** planted-partition and flat/nested ranked models agree
with graph-tool 3.9 to machine precision, synthetic
assortative/non-assortative and ordered/unordered model-selection tests
pass, ranked layouts render cleanly, and both cookbook examples run end
to end. **Complete.**

### 9. Large networks

Build parallel and scalable inference infrastructure for very large
networks (`N > 10^4`, `E > 10^5`): portable header-only C++23
`std::thread` pool (`src/parallel.h`, working identically on macOS Apple
clang without `libomp`, Linux GCC/Clang, and Windows), thread
configuration API
([`rgraphtool_max_threads()`](https://diogro.github.io/rgraphtool/reference/rgraphtool_threads.md),
[`rgraphtool_num_threads()`](https://diogro.github.io/rgraphtool/reference/rgraphtool_threads.md),
[`rgraphtool_set_num_threads()`](https://diogro.github.io/rgraphtool/reference/rgraphtool_threads.md),
`options(rgraphtool.nthreads = ...)`), lock-free read-only parallel
`merge_sweep()` candidate evaluation and parallel batch `mh_sweep()`
screening, thread-safe `LogQ` with larger `MathTables` (`2M` entries),
dynamic-capacity upper `BlockLevel`s with sparse `CompactHierarchy`
snapshots (eliminating `O(L * N)` upper-level overhead), and exact
`O(deg)` incremental upper-level updates for `NestedOverlap`
(`src/overlap.cpp`, matching full recomputation to `< 3.9e-14`).

Implement the guide’s strictly agglomerative search
(`control = list(bisection = FALSE, r = 0.5)`), hierarchy-shape
initialization (`simple_init = TRUE`, `simple_init_reverse`,
`base_factor`, `B_min_base`, `top_factor`, `strict`), and warm-start
support (`init_blocks`, `init_levels`) on
[`minimize_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/minimize_blockmodel_dl.md),
[`refine_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/minimize_blockmodel_dl.md),
and
[`minimize_nested_blockmodel_dl()`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md).

Run and benchmark large-network examples (`demo("large_networks")`)
across 10-core Apple Silicon and multi-core x86_64 Linux (up to 64
threads), documenting hardware, peak memory, and practical scaling
limits.

Benchmarks on modular sparse multigraphs (`<k> = 10`): - **10-core Apple
Silicon (macOS, 32 GB RAM)**: - `N = 10,000, E = 49,794`: flat default
`3.88 s` (1 thr) -\> `0.905 s` (8 thr, `4.3x`, `2.8x` faster than
10-thread graph-tool 3.9’s `2.54 s` with 93 nats lower DL); flat
strictly agglomerative (`r = 0.5`) `0.69 s` (1 thr) -\> `0.249 s` (8
thr); nested default `7.62 s` (1 thr) -\> `1.85 s` (8 thr, `11.1x`
faster than 10-thread graph-tool 3.9’s `20.60 s`); nested `simple_init`
`1.34 s` (1 thr) -\> `0.534 s` (8 thr). - `N = 50,000, E = 249,566`:
flat agglom `7.72 s` (1 thr) -\> `2.14 s` (8 thr, `3.6x`); nested
`simple_init` `16.39 s` (1 thr) -\> `5.54 s` (8 thr, `3.0x`). -
`N = 100,000, E = 499,152`: flat agglom `5.48 s` (8 thr); nested
`simple_init` `14.82 s` (8 thr). Peak memory footprint holding both the
50k and 100k graphs plus `LogQ(5000)` (~100 MB) simultaneously:
`604 MB`. - **Multi-core x86_64 Linux server (GCC 11.4, R 4.5.2, tested
up to 64 threads)**: - `N = 50,000, E = 249,566`: flat agglom `14.53 s`
(1 thr) -\> `6.43 s` (4 thr) -\> `3.46 s` (16 thr, `4.2x`) -\> `3.11 s`
(32 thr) -\> `2.99 s` (64 thr, `4.9x`); nested `simple_init` `31.16 s`
(1 thr) -\> `13.78 s` (4 thr) -\> `8.51 s` (16 thr, `3.7x`). -
`N = 100,000, E = 499,152`: flat agglom `32.38 s` (1 thr) -\> `9.78 s`
(16 thr) -\> `7.90 s` (64 thr, `4.1x`); nested `simple_init` `71.23 s`
(1 thr) -\> `26.70 s` (16 thr) -\> `22.94 s` (64 thr, `3.1x`). Peak RSS
holding both 50k and 100k graphs and 64 thread scratch states:
`791 MB`. - **Practical scaling & limits**: `4`–`16` threads capture
most of the parallel speedup (`3x`–`4.3x`) before sequential
priority-queue merge commits dominate (Amdahl’s law), with incremental
gains up to `32`–`64` threads at `N >= 100,000`. On memory-constrained
machines, lowering `options(rgraphtool.q_cache_max = 1000)` trims ~96 MB
from the partition-count cache.

**Done when:** portable multithreading, compact nested states,
incremental `NestedOverlap` updates, strictly agglomerative search, and
`simple_init` pass correctness and scaling benchmarks up to
`N = 100,000, E = 500,000`, and `demo("large_networks")` runs end to
end. **Complete.**

### 10. Uncertain network reconstruction

Implement native C++ uncertain-network reconstruction engines
(`src/dyn_sbm.h`, `src/uncertain.cpp`, `src/closure.cpp`) with exact
dynamic flat and nested SBM priors supporting \\O(1)\\ / \\O(L)\\ edge
multiplicity edits and \\O(\deg(v))\\ vertex block moves:

- Measured networks with unknown Beta-prior or fixed homogeneous FN/FP
  error rates
  ([`measured_block_state()`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md),
  [`measured_error_posterior()`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md)),
  matching `MeasuredBlockState` (Peixoto 2018).
- Heterogeneous per-dyad Beta-Binomial error models with joint
  hyperparameter MCMC
  ([`mixed_measured_block_state()`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md)),
  matching `MixedMeasuredBlockState` (Peixoto 2018).
- Extraneous edge probability estimates
  ([`uncertain_block_state()`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md)),
  matching `UncertainBlockState` (Peixoto 2018).
- Latent Poisson multigraphs: fast \\O(V + E)\\ EM
  ([`latent_multigraph()`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md))
  and full Bayesian multigraph SBM
  ([`latent_multigraph_block_state()`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md)),
  matching
  [`latent_multigraph()`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md)
  and `LatentMultigraphBlockState` (Peixoto 2020).
- Latent triadic closures
  ([`latent_closure_block_state()`](https://diogro.github.io/rgraphtool/reference/latent_closure_block_state.md),
  [`closure_description_length()`](https://diogro.github.io/rgraphtool/reference/latent_closure_block_state.md),
  [`infer_triadic_closure()`](https://diogro.github.io/rgraphtool/reference/latent_closure_block_state.md)),
  supporting \\L \ge 1\\ generations of zero-inflated triadic closures
  over flat or nested SBM substrates and per-vertex \\(t, m, c)\\
  transitivity statistics, matching `LatentClosureBlockState` (Peixoto
  2022).

Add exact conditional posterior edge probabilities
([`edge_probabilities()`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md)),
joint Bayesian reconstruction
([`reconstruct_network()`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md)),
and the declarative callback-free edge-prediction and evaluation
workflow
([`predict_edges()`](https://diogro.github.io/rgraphtool/reference/predict_edges.md),
[`evaluate_edge_prediction()`](https://diogro.github.io/rgraphtool/reference/predict_edges.md)
with tie-aware ROC-AUC, PR-AUC, Brier score, and ROC/PR curves), plus
runnable cookbook demos (`demo("uncertain_reconstruction")`,
`demo("triadic_closure")`).

Exact description lengths, Beta error posteriors, conditional edge
log-probabilities, and \\L = 2\\ triadic-closure layer entropies and
\\(M, E)\\ triad counts match pinned graph-tool 3.9 to `< 2.4e-12`
across flat and nested, directed and undirected, degree-corrected and
non-degree-corrected configurations; incremental edge, vertex, and
closure moves match full recomputation to `< 3.2e-14`, and the edge MCMC
sampler reproduces the exact enumerated Boltzmann stationary
distribution on small graphs to within Monte Carlo error (`< 0.012`).
Note two deliberate mathematical fixes relative to graph-tool 3.9: (1)
flat SBM edge deltas include exact `edges_dl` changes so incremental
`dS` matches `S(after) - S(before)` identically, and (2)
[`uncertain_block_state()`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md)
uses the exact logit `log(q_default) - log1p(-q_default)` for
`q_default > 0` (fixing upstream’s `log1p(+q_default)` typo). Both
cookbook demos run end to end offline in `< 0.5 s`
(`demo("uncertain_reconstruction")` achieves repeated-trial ROC-AUC
`0.99` with posterior error rates `p = 0.196, q = 0.050` vs true
`0.20, 0.05`, and single-snapshot ROC-AUC `0.91` vs `0.84` raw;
`demo("triadic_closure")` separates closure vs substrate edges with
ROC-AUC `1.00`).

**Done when:** all six uncertain/latent models agree with pinned
graph-tool 3.9 to machine precision, incremental moves and exact
small-graph Boltzmann stationary distributions pass, edge-prediction
evaluation and triadic-closure recovery work end to end, and both
cookbook demos run cleanly. **Complete.**

### 11. Reconstruction from dynamics and behavior

Implement native C++ indirect network reconstruction from dynamics and
behavior (`src/dynamics.cpp`, `R/dynamics.R`), covering all 10
`graph-tool` dynamical and behavioral models with exact \\O(M)\\
incremental local-field updates, multi-cascade support, and flat/nested
SBM priors:

- Ising equilibrium & Glauber dynamics
  ([`pseudo_ising_block_state()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md),
  [`ising_glauber_block_state()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md),
  including 3-state `has_zero = TRUE`) and continuous Ising
  ([`pseudo_cising_block_state()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md),
  [`cising_glauber_block_state()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md)).
- Compartmental epidemic cascades
  ([`epidemics_block_state()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md),
  supporting SI/SIS/SIR/SIRS and `exposed = TRUE` SEI/SEIR, plus
  [`dynamics_transmission_probs()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md)).
- Continuous Gaussian equilibrium & dynamics
  ([`pseudo_normal_block_state()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md)
  with strict diagonal-dominance SPD projection and
  [`dynamics_precision_matrix()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md),
  [`normal_glauber_block_state()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md),
  and
  [`linear_normal_block_state()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md)).
- Generalized \\q\\-state Potts / categorical behavior
  ([`pseudo_potts_block_state()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md),
  [`potts_glauber_block_state()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md)).
- L1 (discretized Laplace) regularization (`xl1`, `tl1`), Gaussian
  weight priors (`normal = TRUE`), and Peixoto (2025 PRX) quantized MDL
  weight and node-parameter priors (`xdist = TRUE`, `tdist = TRUE`),
  plus declarative warm-started L1 regularization paths
  ([`reconstruct_l1_path()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md)),
  greedy MAP
  - finite-temperature posterior MCMC
    ([`reconstruct_dynamics()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md)),
    and native C++ forward simulators
    ([`simulate_ising()`](https://diogro.github.io/rgraphtool/reference/simulate_dynamics.md),
    [`simulate_epidemics()`](https://diogro.github.io/rgraphtool/reference/simulate_dynamics.md),
    [`simulate_normal_dynamics()`](https://diogro.github.io/rgraphtool/reference/simulate_dynamics.md),
    [`simulate_potts()`](https://diogro.github.io/rgraphtool/reference/simulate_dynamics.md)).

Run the synthetic and empirical indirect reconstruction cookbook demos
(`demo("dynamics_reconstruction")`, `demo("empirical_dynamics")`) and
document per-model input and preprocessing requirements (log-returns and
standardization for Gaussian financial time series, compartment codes
for epidemics, spin/category ranges for Ising and Potts).

Exact negative log-likelihoods across all 10 models and all prior
components (geometric edge-count prior, conditioned and unconditioned
L1/Laplace priors, and quantized MDL `xdist`/`tdist` priors) match
pinned graph-tool 3.9 to `< 1.2e-13` (note two intentional behavior
improvements over 3.9:
[`epidemics_block_state()`](https://diogro.github.io/rgraphtool/reference/dynamics_block_state.md)
fixes upstream’s unbound `directed` variable typo in
`EpidemicsBlockState.__init__`, and Potts constructors do not silently
mutate or rescale the user’s `f` interaction matrix). Incremental edge,
node-parameter, and SBM moves match full recomputation to `< 2.1e-13`
across all 10 models. Both cookbook demos execute end to end offline:
`demo("dynamics_reconstruction")` (`7.0 s`) achieves
`ROC-AUC = 1.000, F1 = 1.000` on equilibrium Ising (`26/26`), directed
SIS epidemics (`20/20`, median `t_hat = 0.450` vs true `0.45`), directed
`LinearNormal` (`36/36`), and \\q = 3\\ `PseudoPotts` (`26/26`) while
demonstrating how SBM+MDL eliminates pure L1’s bias-sparsity tradeoff
(`F1 = 1.000` vs `0.26`–`0.42`); `demo("empirical_dynamics")` (`0.31 s`)
reconstructs a strictly SPD partial correlation network on
[`datasets::EuStockMarkets`](https://rdrr.io/r/datasets/EuStockMarkets.html)
(`1,859` trading days) and a 15-asset multi-sector Gaussian network
(`20/20` edges).

**Done when:** all 10 dynamical/behavioral models and L1/quantized
priors agree with pinned graph-tool 3.9 to machine precision, synthetic
and empirical reconstruction workflows and L1 regularization paths run
end to end, and input preprocessing is documented. **Complete.**

### 12. Animation and interactive visualization

Provide R-native equivalents for interactive animation, SIRS epidemics,
dynamic layouts, and interactive visualizations (`src/layout.cpp`,
`R/animation.R`), running simulation and layout loops natively in C++
with zero external R package, GTK+, or Cairo dependencies:

- Native Hu (2005) / graph-tool SFDP force-directed layout
  ([`sfdp_layout()`](https://diogro.github.io/rgraphtool/reference/sfdp_layout.md))
  with Barnes-Hut quadtree acceleration (\\N \> 128\\), warm-starting
  (`pos`, `max_iter = 1L`), vertex pinning (`pin`), and
  inter-/intra-group forces (`groups`, `gamma`, `mu`).
- Exact directed/undirected breadth-first search traversal and
  spanning-tree extraction
  ([`graph_bfs()`](https://diogro.github.io/rgraphtool/reference/sfdp_layout.md)).
- Base-R graph rendering
  ([`plot_graph()`](https://diogro.github.io/rgraphtool/reference/plot_graph.md))
  supporting vertex/edge visibility filtering, highlight halos, directed
  arrowheads, labels, and vector/raster R graphics devices.
- Declarative SIRS epidemic animation
  ([`animate_sirs()`](https://diogro.github.io/rgraphtool/reference/animate_sirs.md),
  matching `animation_sirs.py` and `animation_zombies.py` via
  `hide_recovered = TRUE`, newly-infected halos, and
  `fear_state = TRUE`), dynamic spatial-segregation layout animation
  ([`animate_dynamic_layout()`](https://diogro.github.io/rgraphtool/reference/animate_sirs.md),
  matching `animation_dancing.py` with incremental SFDP +
  distance-shortening rewiring), and interactive BFS shortest-path tree
  explorer
  ([`interactive_bfs()`](https://diogro.github.io/rgraphtool/reference/interactive_bfs.md),
  matching `interactive_bst.py`).

Run all three cookbook demos end to end offline
(`demo("animation_sirs")`, `demo("animation_dynamic_layout")`,
`demo("interactive_bfs")`) and document platform and display
requirements in Rd help pages:

- `mode = "html"` generates zero-dependency standalone HTML5 Canvas
  viewers that run headless during creation and interactively offline in
  any modern browser or RStudio Viewer (play/pause, frame scrubber, FPS
  control, live sparklines, and hover/click BFS inspection).
- `mode = "frames"` exports numbered PNG frames via
  [`grDevices::png()`](https://rdrr.io/r/grDevices/png.html) (standard
  macOS/Windows/Cairo-enabled Linux R builds) for external video or GIF
  encoding.
- `mode = "device"` and `mode = "locator"` require an interactive R
  graphics device (Quartz, X11, Windows, or RStudioGD); `locator` errors
  cleanly in non-interactive sessions.

[`graph_bfs()`](https://diogro.github.io/rgraphtool/reference/sfdp_layout.md)
matches pinned graph-tool 3.9 `shortest_distance(..., max_dist=3)`
identically across directed and undirected test graphs, and all 78
animation and layout unit tests verify determinism under
[`set.seed()`](https://rdrr.io/r/base/Random.html), empty/singleton/
self-loop/parallel-edge/disconnected handling, quadtree execution,
HTML5/PNG output generation, and spatial-segregation edge-length
contraction. All three cookbook demos complete offline in `< 0.3 s`
combined (`demo("animation_sirs")` in `0.13 s`,
`demo("animation_dynamic_layout")` in `0.09 s` with `25.7%` mean
edge-length contraction over 80 steps, and `demo("interactive_bfs")` in
`0.04 s`).

**Done when:** SFDP layout, BFS traversal, SIRS/fear epidemic
animations, spatial-segregation dynamic layouts, and interactive BFS
visualization work end to end offline with documented platform and
display requirements. **Complete.**

### 13. Plotting integration

Provide R plotting and map-integration equivalents for the matplotlib
and basemap examples (`R/plot-integration.R`, `R/plot-map.R`, enhanced
[`plot_graph()`](https://diogro.github.io/rgraphtool/reference/plot_graph.md)),
supporting both vector (`pdf`, `svg`, `ps`, `eps`) and raster (`png`,
`jpeg`, `tiff`) output without requiring Python:

- Property-to-size scaling
  ([`prop_to_size()`](https://diogro.github.io/rgraphtool/reference/plotting_integration.md)),
  HCL colour mapping and standalone colorbars
  ([`map_colors()`](https://diogro.github.io/rgraphtool/reference/plotting_integration.md),
  [`plot_colorbar()`](https://diogro.github.io/rgraphtool/reference/plotting_integration.md)),
  curved edges, vertex/edge z-ordering, axes/coordinate grids, and
  unified multi-format export
  ([`save_graph_plot()`](https://diogro.github.io/rgraphtool/reference/plotting_integration.md),
  using `svglite` when installed for SVG without XQuartz/Cairo
  dependencies and falling back to
  [`grDevices::svg()`](https://rdrr.io/r/grDevices/cairo.html)).
- Embeddable base-R `grid` grobs
  ([`graph_grob()`](https://diogro.github.io/rgraphtool/reference/plotting_integration.md)),
  tidy vertex/edge data frames
  ([`graph_to_data_frames()`](https://diogro.github.io/rgraphtool/reference/plotting_integration.md)),
  and optional `ggplot2` integration
  ([`ggplot_graph()`](https://diogro.github.io/rgraphtool/reference/plotting_integration.md)).
- Pure-R cartographic projections
  ([`project_coords()`](https://diogro.github.io/rgraphtool/reference/plot_graph_map.md):
  Robinson, Mercator, equirectangular, Mollweide, and orthographic
  hemisphere clipping), spherical great-circle geodesic edges
  ([`great_circle_arc()`](https://diogro.github.io/rgraphtool/reference/plot_graph_map.md),
  splitting cleanly across the antimeridian), built-in offline world
  coastlines
  ([`world_coastlines()`](https://diogro.github.io/rgraphtool/reference/plot_graph_map.md),
  plus custom polygon support), and geographic network plotting
  ([`plot_graph_map()`](https://diogro.github.io/rgraphtool/reference/plot_graph_map.md)).

Run both cookbook demos end to end offline without Python
(`demo("plotting_integration")` in `0.91 s` exporting multi-panel
PDF/PNG, vector SVG, `grid`-composed PDF, and `ggplot2` PNG;
`demo("plotting_map")` in `0.31 s` exporting Robinson PDF/PNG,
Orthographic globe PNG, and regional Mercator PDF).

All 75 plotting-integration and map tests pass (`0` failures, `0`
warnings), verifying exact
[`prop_to_size()`](https://diogro.github.io/rgraphtool/reference/plotting_integration.md)
scaling, Mollweide and orthographic projection geometry, trans-Pacific
antimeridian great-circle splitting, `grid` grob and `ggplot2`
rendering, and binary file headers across PDF, SVG, EPS, PNG, JPEG, and
TIFF outputs.

**Done when:** base-R multi-panel, `grid`, optional `ggplot2`, and
geographic map-integration workflows run end to end offline across
vector and raster formats without Python. **Complete.**

### 14. C++ extensions

Document and test a public R/C++ extension interface
(`inst/include/rgraphtool.h`, moving vendored `boost::adj_list` headers
to `inst/include/rgraphtool/vendor/` so `src/graph_types.h` and external
C++ extensions share a single ABI definition via
`LinkingTo: rgraphtool, Rcpp, BH` or
[`Rcpp::sourceCpp()`](https://rdrr.io/pkg/Rcpp/man/sourceCpp.html)):

- Zero-copy graph extraction (`rgraphtool::get_graph(SEXP)` accepting
  either an `rgraphtool_graph` R object or its external pointer) and
  inspection helpers (`num_vertices`, `num_edges`, `is_directed`,
  `source`, `target`, `edge_index`, `out_degree`, `in_degree`,
  `total_degree`, `checked_vertex`).
- Range-based graph iteration (`vertices_range`, `edges_range`,
  `out_edges_range`, `in_edges_range`, `all_edges_range`,
  `out_neighbors_range`, `in_neighbors_range`, `all_neighbors_range`).
- Vertex and edge property maps (`vertex_property_view`,
  `edge_property_view`, `make_vertex_property`, `make_edge_property`)
  with per-access bounds-checked defaults (`VertexPropertyMap<T>`,
  `EdgePropertyMap<T>`, `ConstEdgePropertyMap<T>`) and upfront-validated
  zero-overhead unchecked views via `.get_unchecked()`
  (`UncheckedVertexPropertyMap<T>`, `UncheckedEdgePropertyMap<T>`,
  `UncheckedConstEdgePropertyMap<T>`).

Provide native R and standalone C++ extension equivalents of the
cookbook examples (`src/extensions.cpp`, `R/extensions.R`,
`man/cpp_extensions.Rd`):

- Built-in
  [`kcore_decomposition()`](https://diogro.github.io/rgraphtool/reference/cpp_extensions.md)
  (exact \\O(N + E)\\ Batagelj-Zaversnik algorithm matching pinned
  graph-tool 3.9
  [`kcore_decomposition()`](https://diogro.github.io/rgraphtool/reference/cpp_extensions.md)
  identically on directed/undirected graphs, self-loops, parallel edges,
  and isolated vertices) and
  [`vertex_strength()`](https://diogro.github.io/rgraphtool/reference/cpp_extensions.md),
  each supporting `checked = TRUE/FALSE`.
- Extension discovery and compilation helpers
  ([`rgraphtool_include_dir()`](https://diogro.github.io/rgraphtool/reference/cpp_extensions.md),
  [`rgraphtool_extension_example()`](https://diogro.github.io/rgraphtool/reference/cpp_extensions.md),
  [`compile_cpp_extension()`](https://diogro.github.io/rgraphtool/reference/cpp_extensions.md))
  plus bundled standalone `.cpp` examples
  (`inst/examples/extensions/kcore_extension.cpp`,
  `inst/examples/extensions/weighted_degree_extension.cpp`) and
  end-to-end demo (`demo("cpp_extensions")`).

All 34 C++ extension tests (`tests/testthat/test-extensions.R`) pass
(`0` failures, `0` warnings), including live compilation and execution
of both bundled `.cpp` files and inline C++ code via
[`compile_cpp_extension()`](https://diogro.github.io/rgraphtool/reference/cpp_extensions.md),
exact parity against pinned graph-tool 3.9 fixtures, and explicit
verification of `std::out_of_range` bounds checks on checked views and
upfront capacity checks on `.get_unchecked()`. `demo("cpp_extensions")`
compiles and executes all three extensions offline in `~6.9 s`,
measuring `0.005 s` unchecked vs. `0.013 s` checked (`2.6x` faster)
across 50 \\k\\-core passes on a 2,000-vertex, 9,996-edge weighted
network.

**Done when:** the public `inst/include/rgraphtool.h` interface is
documented and tested for graph iteration, property map access, and
checked vs. unchecked access, and standalone C++ extension examples
compile and run end to end. **Complete.**

## Completion requirements for every milestone

- Add runnable, documented R examples corresponding to the relevant
  cookbook sections, with data sources, attribution, and reproducible
  setup.
- Verify numerical behavior against a pinned official graph-tool
  reference, using hand-checkable cases and regression fixtures as
  appropriate. Record deliberate differences and stochastic tolerances.
- Test invalid inputs, graph direction, empty graphs, loops, and
  parallel edges where applicable; run the relevant package tests and
  `R CMD check`.
- Specify runtime and memory limits, random-seed behavior, and
  platform-specific dependencies. Do not silently substitute a different
  statistical objective.
- Mark a milestone complete only when its examples execute end to end; a
  partial API or scaffold is not cookbook parity.
