# Configure parallel worker threads for large-network inference

Query and configure the portable C++23 thread pool used to parallelize
large-network stochastic block model minimization in
[`minimize_blockmodel_dl`](https://diogro.github.io/rgraphtool/reference/minimize_blockmodel_dl.md),
[`refine_blockmodel_dl`](https://diogro.github.io/rgraphtool/reference/minimize_blockmodel_dl.md),
and
[`minimize_nested_blockmodel_dl`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md).

## Usage

``` r
rgraphtool_max_threads()
rgraphtool_num_threads()
rgraphtool_set_num_threads(nthreads)
```

## Arguments

- nthreads:

  Positive integer number of worker threads (capped at 256), or `"auto"`
  to use `rgraphtool_max_threads()`.

## Details

rgraphtool uses a built-in C++23 `std::thread` pool rather than OpenMP,
so multithreading works identically on macOS (including Apple clang
without `libomp`), Linux, and Windows with no external threading
library. Worker threads never call the R API or R's RNG; per-thread
deterministic `xoshiro256**` generators are seeded from R's RNG on the
main thread before entering parallel regions, so results remain
reproducible under [`set.seed()`](https://rdrr.io/r/base/Random.html)
for a fixed thread count.

By default rgraphtool runs on `1L` thread so background tasks and shared
systems are never oversubscribed unexpectedly. Opt in globally via
`rgraphtool_set_num_threads(k)` or `options(rgraphtool.nthreads = k)`,
or per call via the `nthreads` argument.

In practice, 4–16 threads typically provide the best speedups (roughly
3x–4.3x on networks with \\N = 10^4\\–\\10^5\\ vertices), after which
sequential priority-queue merge commits limit further scaling; 32–64
threads still give incremental gains on \\N \ge 10^5\\ networks. To
reduce peak RAM on very large networks, lower
`options(rgraphtool.q_cache_max = 1000)` (the default `5000` uses ~100
MB for the exact partition-count table) and pair multithreading with
strictly agglomerative search
(`control = list(bisection = FALSE, r = 0.5)`) or hierarchy-shape
initialization (`simple_init = TRUE`).

## Value

`rgraphtool_max_threads()` returns the number of hardware threads
reported by `std::thread::hardware_concurrency()` (at least `1L`).

`rgraphtool_num_threads()` returns the active default thread count.

`rgraphtool_set_num_threads()` invisibly returns the configured integer
thread count.

## See also

[`minimize_blockmodel_dl`](https://diogro.github.io/rgraphtool/reference/minimize_blockmodel_dl.md),
[`minimize_nested_blockmodel_dl`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md),
`demo("large_networks")`

## Examples

``` r
rgraphtool_max_threads()
#> [1] 10
old <- rgraphtool_num_threads()
rgraphtool_set_num_threads(min(2L, rgraphtool_max_threads()))
rgraphtool_num_threads()
#> [1] 2
rgraphtool_set_num_threads(old)
```
