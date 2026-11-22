# graph-tool-R roadmap

## Goal and scope

Build a native R/C++ port that can run R equivalents of the
[graph-tool cookbook](https://graph-tool.skewed.de/static/docs/stable/demos/index.html),
following its order from simpler to more advanced examples. Start with block
models and the inference guide. Python is not a runtime dependency.

The cookbook targets current graph-tool. Use the
[official upstream](https://git.skewed.de/count0/graph-tool) and a pinned
release or commit as the reference for each milestone; record any version
differences. The existing adjacency headers came from a 2020 checkout and
must not be treated as a current inference backend. Preserve upstream license
notices and attribution when adapting code or distributing data.

R examples should express equivalent functionality, rather than copying
Python syntax. Clearly document differences in indexing, mutable state,
defaults, output formats, and supported graph types. The matplotlib and C++
extension chapters need R-native equivalents, not a Python compatibility
layer.

## Current state

- [x] Native directed and undirected multigraph construction and inspection.
- [x] Fixed-partition block sizes, edge-count matrix, and conditional
  microcanonical likelihood.
- [x] Flat block-model description length and experimental greedy fitting.
- [ ] Hierarchical inference, posterior sampling, cookbook datasets, and
  visualization. The greedy search is not the upstream multilevel optimizer.

Fixed-partition scores have been checked against graph-tool 2.92, including a
football-network partition. In one seeded football run, the greedy fit scored
1869.62 nats versus 1795.95 for an upstream fit; these stochastic results
illustrate the current optimizer's limitations, not a parity guarantee.

## Milestones

### 1. First flat SBM example

- [x] Port and validate the flat microcanonical SBM description length,
  including partition, edge-count, and degree-sequence priors. Support both
  degree-corrected and non-degree-corrected models, directed and undirected
  graphs, parallel edges, and self-loops.
- [ ] Implement scalable native description-length minimization and expose fitted block
  memberships, block sizes, edge counts, and description length. Document the
  search algorithm and its differences from upstream; do not promise a global
  optimum. An initial greedy vertex-move and merge search is available, but
  is slower and can find worse optima than the upstream multilevel optimizer.
- [ ] Provide a reproducible way to obtain and load the cookbook's football
  graph, including its layout or an explicitly documented R equivalent.
- [ ] Provide an R example that fits the football graph, inspects memberships
  and block counts, and produces the partition and affinity-matrix plots.

**Done when:** the example runs end to end; fixed-partition scores agree with
the pinned upstream reference within stated numerical tolerances; tests cover
edge cases and compare inferred objective values without requiring identical
stochastic labels or partitions.

### 2. Hierarchical partitions and refinement

- [ ] Implement nested block states, hierarchical description length,
  minimization, level traversal, and hierarchy summaries.
- [ ] Load the C. elegans neural network and produce the hierarchy visualization
  and per-level membership examples.
- [ ] Add merge-split refinement and annealing, with reproducible random seeds.

**Done when:** the hierarchical cookbook examples run and refinement tests
verify that zero-temperature accepted moves do not worsen the objective.

### 3. Model selection

- [ ] Support comparing degree-corrected and non-degree-corrected nested models,
  including description-length differences and log posterior odds.
- [ ] Run the political blogs and football model-selection examples.

**Done when:** scores agree with the reference for identical states and the R
examples report their comparisons without assuming stochastic fits are identical.

### 4. Posterior sampling

- [ ] Implement the flat and nested MCMC moves and equilibration controls used
  by the guide, including merge-split sampling.
- [ ] Add partition collection, label alignment, vertex marginals, consensus
  estimates, and posterior mode characterization.
- [ ] Run the Les Misérables and hierarchical sampling examples with the
  corresponding visualizations.

**Done when:** exact small-graph checks validate stationary probabilities or
transition behavior, and seeded integration examples run with sampling
uncertainty and convergence limitations documented.

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
