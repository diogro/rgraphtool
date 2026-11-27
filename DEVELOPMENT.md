# Developing `rgraphtool`

Guide for compiling, testing, benchmarking against upstream `graph-tool`, and
extending `rgraphtool`.

## 1. Prerequisites & Build Toolchain

- **R**: `R >= 4.1.0`.
- **C++ compiler**: C++23-capable compiler (`g++` / `clang++` supporting
  `-std=c++23` or `CXX23`), needed to compile `src/*.cpp`.
- **Required R packages**: `Rcpp`, `BH (>= 1.81.0)`, `grid`.
- **Suggested R packages**: `testthat (>= 3.0.0)`, `igraph`, `ggplot2`,
  `svglite`, `knitr`, `rmarkdown`, `roxygen2`, `pkgdown`.
- **Upstream reference parity (optional)**: `pixi` environment in
  `tools/reference/` running upstream `conda-forge` `graph-tool` 3.9 for
  entropy and partition parity checks.
- Set `R_LIBS_USER` if using a custom user library:

```bash
export R_LIBS_USER=/dados/home/melo/R/x86_64-pc-linux-gnu-library/4.5
```

## 2. Repository Layout

```text
rgraphtool/
├── src/                # C++23 cores, Rcpp bindings, Makevars(.win), parallel.h
├── inst/include/       # Public C++ headers exposed to LinkingTo extensions
├── R/                  # R S3 wrappers, validation, plotting, data loaders
├── man/                # Generated Rd files (via roxygen2)
├── tests/testthat/     # Unit, parity, dynamics, and performance tests
├── vignettes/          # R-native graph-tool cookbook chapters
├── demo/               # Interactive/benchmark demos
├── tools/reference/    # Upstream graph-tool 3.9 pixi verification helpers
├── pkgdown/            # Website build script and theme CSS
├── AGENTS.md           # Contributor/agent invariants & C++ pitfalls
└── DEVELOPMENT.md      # This file
```

## 3. Build, Test, Check, and Website Workflows

### Regenerate Rcpp wrappers & Rd docs

```bash
Rscript -e 'Rcpp::compileAttributes("."); roxygen2::roxygenise(".")'
```

### Compile and install locally

```bash
# Multi-core C++ compile via MAKEFLAGS
MAKEFLAGS="-j$(nproc)" R CMD INSTALL --library="$R_LIBS_USER" .
```

### Run tests

```bash
# Run a specific context (e.g. topology, planted-partition, dynamics)
Rscript -e 'library(rgraphtool); testthat::test_dir("tests/testthat", filter = "^topology$", stop_on_failure = TRUE)'

# Full test suite
Rscript -e 'library(rgraphtool); testthat::test_dir("tests/testthat", stop_on_failure = TRUE)'
```

### Source tarball & `R CMD check`

Run tarball build and check outside the git root in `/tmp`:

```bash
mkdir -p /tmp/rgt_check && (
  cd /tmp/rgt_check &&
  R CMD build /home/melo/projects/rgraphtool &&
  R CMD check --no-manual rgraphtool_*.tar.gz
) && rm -rf /tmp/rgt_check
```

Target: `Status: OK` (`0 errors`, `0 warnings`, `0 notes`).

### Rebuild the `pkgdown` site

```bash
Rscript pkgdown/build_site.R         # full rebuild into docs/
Rscript pkgdown/build_site.R --lazy  # incremental rebuild
Rscript -e 'pkgdown::check_pkgdown()'
```

## 4. Extending `rgraphtool` Safely

1. **Validate in R, execute in C++**: Keep argument checks (`validate_count`,
   vertex bounds, finite weights, matrix dimensions) in `R/` with clear errors
   (`call. = FALSE`), and implement heavy loops/sampling in C++23 inside `src/`.
2. **Thread safety**: Spawn worker threads via `src/parallel.h` respecting
   `rgraphtool_num_threads()`, keep R's `Rcpp` / `SEXP` objects and RNG isolated
   on the main thread, and pass plain C++ containers/views into worker pools.
3. **Public C++ API**: If exposing new internal graph/state structs to
   extensions, update `inst/include/` headers and verify
   `compile_cpp_extension()` still builds cleanly.
4. **Parity tests**: When modifying SBM entropy (`blockmodel`, `weighted`,
   `layered`, `planted_partition`, `ranked`, `bipartite`, `dynamics`), add
   exact numerical parity assertions comparing against closed-form identities or
   upstream `graph-tool` 3.9 reference values.

## 5. Pre-Commit Checklist

- [ ] `Rcpp::compileAttributes()` and `roxygen2::roxygenise()` are up to date.
- [ ] `testthat::test_dir("tests/testthat", stop_on_failure = TRUE)` passes with
      `0` failures and `0` warnings.
- [ ] `R CMD check --no-manual` reports `Status: OK`.
- [ ] `pkgdown/build_site.R` and `pkgdown::check_pkgdown()` succeed cleanly.
- [ ] Temporary build artifacts (`.o`, `.so` in untracked dirs if building clean,
      `/tmp/rgt_check`) are clean and root helper docs remain in `.Rbuildignore`.
