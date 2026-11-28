# Agent & Contributor Guide (`rgraphtool`)

Rules, C++23 architecture notes, and upstream `graph-tool` parity invariants
for working safely inside `rgraphtool`.

## 1. Repository Boundary & Upstream Target

- **Role**: Domain-agnostic native R/C++23 port of Python `graph-tool` release
  **3.9** (upstream commit `b24bc428e8426d0acb36b2d57bd61633335a3a26`), with zero
  Python runtime dependencies.
- **Keep domain logic out**: Transcriptomics workflows belong in downstream
  packages such as `gensbm`. Keep `rgraphtool` focused on general graphs, SBM
  inference/MCMC, dynamics/network reconstruction, topology/centrality,
  layouts, visualization, and C++ extension headers (`inst/include/`).
- **License provenance**: Licensed under LGPL (>= 3). Preserve upstream
  `graph-tool` / Boost copyright and license headers whenever porting or
  extending C++ headers/sources.
- **Do not commit or push automatically** unless explicitly requested by the user.

## 2. Module Map (`src/` and `R/`)

- `src/graph.cpp`, `src/graph_types.h`, `R/graph.R`, `R/io.R`: Core directed /
  undirected multigraph storage, Boost adjacency wrappers, vertex/edge/graph
  properties, `.gt` binary I/O, and `igraph` conversion.
- `src/parallel.h`, `R/parallel.R`: Portable C++23 `std::thread` pool and
  thread configuration (`rgraphtool_set_num_threads()`).
- `src/blockmodel.cpp`, `src/weighted.cpp`, `src/weight_cov.h`, `src/layered.cpp`,
  `src/planted_partition.cpp`, `src/ranked.cpp`, `src/bipartite.cpp`,
  `src/overlap.cpp`: Exact MDL entropy, conditional likelihoods, and
  agglomerative minimization / merge-split sweeps (`R/block-state.R`,
  `R/nested.R`, `R/weighted.R`, `R/layered.R`, `R/planted-partition.R`,
  `R/ranked.R`, `R/bipartite.R`, `R/overlap.R`, `R/mcmc.R`,
  `R/model-selection.R`, `R/partition-modes.R`).
- `src/dynamics.cpp`, `src/dyn_sbm.h`, `src/closure.cpp`, `src/uncertain.cpp`:
  Network reconstruction from Ising, Potts, Glauber, epidemic, and Gaussian /
  PseudoNormal dynamics (`R/dynamics.R`, `R/closure.R`, `R/uncertain.R`).
- `src/topology.cpp`, `R/topology.R`: Centrality, shortest paths, components,
  clustering, assortativity, modularity/Louvain, and random graph generators.
- `src/layout.cpp`, `R/plot-*.R`, `R/animation.R`: SFDP layouts, hierarchy /
  bipartite / marginal / map plotting, and HTML5 Canvas animations.
- `inst/include/` & `src/extensions.cpp`: Public C++ headers for downstream
  `LinkingTo: rgraphtool` packages (`compile_cpp_extension()`).

## 3. Critical Pitfalls & Invariants Checklist

1. **Edge order inside Boost adjacency lists**:
   `graph_edges()` is source-major: edges are grouped by source vertex, keep
   insertion order within a source (not sorted by target), and keep the
   orientation they were given. `graph_add_edges()` reorders attached edge
   properties to match, and positions change after every add.
   - Any external weight vector passed to `assortativity()`, `modularity()`, or
     block states **must** align with the current `graph_edges(g)`, or else be
     attached via `graph_add_edges(..., properties = list(weight = ...))` /
     `graph_set_edge_property()` and referenced by property name (`weights = "weight"`).
   - Numeric edge properties live in C++, all others in the R environment; a
     property may not change between the two (it is an error).
2. **1-based R vs. 0-based C++ vertex indices**:
   R exposes 1-based vertex IDs (`1..N`) and block labels (`1..B`), while C++
   uses 0-based indices (`it->s`, `it->t`). Validate lengths, finiteness, and
   `1..N` bounds in R wrappers before crossing `.Call()` / `Rcpp`.
3. **Regenerate `RcppExports` whenever modifying `[[Rcpp::export]]`**:
   If you add, rename, or change the C++ signature of any `// [[Rcpp::export]]`
   function in `src/*.cpp`, always run:
   ```r
   Rcpp::compileAttributes(".")
   ```
   before `R CMD INSTALL` or `R CMD check`. Do not hand-edit `src/RcppExports.cpp`
   or `R/RcppExports.R`. The Rd files in `man/` and `NAMESPACE` are maintained
   **by hand** (not by roxygen): after any change to an exported signature,
   update the Rd `\usage` / `\arguments` and run `tools::codoc(dir = ".")`.
4. **One log q(n, k) table**: every SBM family takes degree-DL partition
   counts from `src/logq.h` (`rgt_logq::shared()`). Never add a model-local
   copy; per-model copies with a 2000 cap and a wrong approximation put
   weighted / layered / PP degree DLs hundreds of nats off graph-tool.
5. **No BLAS on pool workers**: BLAS/LAPACK may be called from the main thread
   only. Inside `rgt_parallel::parallel_for` use plain C++ kernels (a
   multithreaded OpenBLAS would oversubscribe cores; some builds are not
   reentrant).
6. **C++23 & portable threading rules**:
   `src/Makevars` / `src/Makevars.win` enforce `CXX_STD = CXX23` and
   `PKG_CPPFLAGS = -I../inst/include`. Never rely on compiler-specific OpenMP
   pragmas for core inference; use `src/parallel.h` (`std::thread`) and avoid
   GNU extensions that fail strict CRAN / portable C++23 checks. Keep R API
   calls (memory allocation, RNG if using R's RNG, exceptions) off worker
   threads.
7. **Undirected half-edge (stub) & self-loop conventions**:
   `blockmodel.cpp`, `planted_partition.cpp`, and `topology.cpp` follow
   graph-tool's convention where undirected internal edges count **2 stubs** in
   $e_{rr}$ (`W = 2E`) and undirected modularity sums $e_{rr}/(2m) - (k_r/(2m))^2$.
   Verify exact parity against Python `graph-tool` 3.9 using the pixi
   environment under `tools/reference/` whenever changing entropy, Log-BF, or
   MCMC acceptance formulas.
8. **Declarative API over closure callbacks**:
   Never require user-supplied R closures inside tight MCMC inner loops. Expose
   declarative parameters (`nsamples`, `collect`, `history`, `refine`,
   `mode_cluster_state()`) executed natively in C++.
9. **Build, check, and site cleanups**:
   - Register exported R functions and datasets in `_pkgdown.yml` when updating
     documentation.
   - Run `R CMD build` / `R CMD check --no-manual` in `/tmp` and `cd` out before
     removing the temporary directory.
   - Keep `AGENTS.md`, `DEVELOPMENT.md`, `ROADMAP.md`, `tools/`, and `docs/`
     listed in `.Rbuildignore`.
