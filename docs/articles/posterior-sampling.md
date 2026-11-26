# Posterior MCMC Sampling & Partition Modes

This article corresponds to the `posterior_sampling` and
`hierarchical_sampling` demos. Rather than relying on a single point
estimate of the community partition, `rgraphtool` samples directly from
the exact Bayesian posterior \\(P(b \mid A) \propto \exp(-\Sigma(A,
b))\\) using merge-split MCMC and summarizes multimodal posteriors
through optimal label alignment and mode clustering.

## 1. Flat posterior sampling, vertex marginals, and mode clustering

Starting from a multilevel point estimate,
[`mcmc_equilibrate()`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md)
runs merge-split MCMC until the description length reaches stationarity,
after which
[`sample_posterior()`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md)
draws thinned posterior partitions.

\
`lesmis`` ``<-`` ``load_or_fallback_lesmis``(``)`\
`g`` ``<-`` ``lesmis``$``graph`\
\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``42``)`\
`state`` ``<-`` `[`minimize_blockmodel_dl`](https://diogro.github.io/rgraphtool/reference/minimize_blockmodel_dl.md)`(``g``)`\
`eq`` ``<-`` `[`mcmc_equilibrate`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md)`(``state``, wait ``=`` ``100L``, nbreaks ``=`` ``2L``,`\
`                       mcmc_args ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``niter ``=`` ``10L``)``)`\
[`cat`](https://rdrr.io/r/base/cat.html)`(`[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"Equilibrated flat DC SBM in %d steps (DL = %.2f nats)\n"``,`\
`            ``eq``$``niter``, ``eq``$``entropy``)``)`\
`#> Equilibrated flat DC SBM in 496 steps (DL = 716.56 nats)`\
\
`post`` ``<-`` `[`sample_posterior`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md)`(``eq``$``state``, nsamples ``=`` ``300L``, thin ``=`` ``8L``, burnin ``=`` ``0L``)`\
\
`# Posterior distribution over the number of groups B`\
`B_samples`` ``<-`` `[`vapply`](https://rdrr.io/r/base/lapply.html)`(``post``$``samples``, ``function``(``b``)`` `[`max`](https://rdrr.io/r/base/Extremes.html)`(``b``)``, `[`integer`](https://rdrr.io/r/base/integer.html)`(``1L``)``)`\
`B_dist`` ``<-`` `[`prop.table`](https://rdrr.io/r/base/proportions.html)`(`[`table`](https://rdrr.io/r/base/table.html)`(``B_samples``)``)`\
[`round`](https://rdrr.io/r/base/Round.html)`(``B_dist``, ``3``)`\
`#> B_samples`\
`#>     5     6     7     8     9 `\
`#> 0.013 0.067 0.723 0.190 0.007`\
\
`# Single-mode optimal label alignment and multimodal clustering`\
`pm`` ``<-`` `[`partition_mode_state`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md)`(``post``)`\
[`print`](https://rdrr.io/r/base/print.html)`(``pm``)`\
`#> <partition mode: M = 300, B = 11, uncertainty = 0.1581, DL = 14761.6343 nats>`\
\
`mc`` ``<-`` `[`mode_cluster_state`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md)`(``post``)`\
[`print`](https://rdrr.io/r/base/print.html)`(``mc``)`\
`#> <flat mode cluster state: K = 3 modes (Be = 1.98), M = 300 partitions, DL = 12859.1903 nats>`\
`#>  mode   M      weight  B uncertainty posterior_entropy description_length`\
`#>     1 194 0.646666667  7   0.1049003          27.76036          6993.0107`\
`#>     2 104 0.346666667  8   0.1727023          35.92085          5321.8706`\
`#>     3   2 0.006666667 10   0.1363636          14.55609           327.4646`\
\
`old_par`` ``<-`` ``graphics``::`[`par`](https://rdrr.io/r/graphics/par.html)`(``mfrow ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``1L``, ``2L``)``, mar ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``4``, ``4``, ``3``, ``1``)``)`\
`graphics``::`[`barplot`](https://rdrr.io/r/graphics/barplot.html)`(``B_dist``, xlab ``=`` ``"Number of groups B"``, ylab ``=`` ``"Posterior P(B | A)"``,`\
`                  col ``=`` ``"steelblue"``, main ``=`` ``"Posterior P(B | A)"``)`\
\
`dom`` ``<-`` ``mc``$``modes``[[``1L``]``]`\
[`plot_marginals`](https://diogro.github.io/rgraphtool/reference/plot_marginals.md)`(``g``, ``dom``, positions ``=`` ``lesmis``$``positions``, radius ``=`` ``0.045``,`\
`               main ``=`` `[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"Dominant mode (w = %.2f, unc = %.3f)"``,`\
`                              ``mc``$``weights``[``1L``]``, `[`mode_uncertainty`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md)`(``dom``)``)``)`

![](posterior-sampling_files/figure-html/flat-posterior-1.png)

\
`graphics``::`[`par`](https://rdrr.io/r/graphics/par.html)`(``old_par``)`

## 2. Hierarchical posterior sampling and consensus hierarchies

The exact same workflow extends to nested SBMs, capturing uncertainty in
both level-0 memberships and the depth and branching structure of the
hierarchy.

\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``42``)`\
`nstate`` ``<-`` `[`minimize_nested_blockmodel_dl`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md)`(``g``)`\
[`print`](https://rdrr.io/r/base/print.html)`(``nstate``)`\
`#> <nested block state: degree-corrected, 8 levels, description length 708.3461 nats>`\
`#> l: 0, N: 77, B: 7`\
`#> l: 1, N: 7, B: 1`\
`#> l: 2, N: 1, B: 1`\
\
`npost`` ``<-`` `[`sample_posterior`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md)`(``nstate``, nsamples ``=`` ``250L``, thin ``=`` ``5L``,`\
`                          equilibrate_args ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``wait ``=`` ``100L``, nbreaks ``=`` ``2L``)``)`\
\
`B0`` ``<-`` `[`vapply`](https://rdrr.io/r/base/lapply.html)`(``npost``$``samples``, ``function``(``bs``)`` `[`max`](https://rdrr.io/r/base/Extremes.html)`(``bs``[[``1L``]``]``)``, `[`integer`](https://rdrr.io/r/base/integer.html)`(``1L``)``)`\
`depth`` ``<-`` `[`vapply`](https://rdrr.io/r/base/lapply.html)`(``npost``$``samples``, ``function``(``bs``)`` ``{`\
`  ``top`` ``<-`` `[`match`](https://rdrr.io/r/base/match.html)`(``1L``, `[`vapply`](https://rdrr.io/r/base/lapply.html)`(``bs``, ``max``, `[`integer`](https://rdrr.io/r/base/integer.html)`(``1L``)``)``)`\
`  ``if`` ``(`[`is.na`](https://rdrr.io/r/base/NA.html)`(``top``)``)`` `[`length`](https://rdrr.io/r/base/length.html)`(``bs``)`` ``else`` ``top`\
`}``, `[`integer`](https://rdrr.io/r/base/integer.html)`(``1L``)``)`\
\
[`round`](https://rdrr.io/r/base/Round.html)`(`[`prop.table`](https://rdrr.io/r/base/proportions.html)`(`[`table`](https://rdrr.io/r/base/table.html)`(``B0``)``)``, ``3``)`\
`#> B0`\
`#>     6     7     8     9 `\
`#> 0.016 0.428 0.352 0.204`\
[`round`](https://rdrr.io/r/base/Round.html)`(`[`prop.table`](https://rdrr.io/r/base/proportions.html)`(`[`table`](https://rdrr.io/r/base/table.html)`(``depth``)``)``, ``3``)`\
`#> depth`\
`#>     2     3     4     5     6     7 `\
`#> 0.036 0.668 0.200 0.068 0.024 0.004`\
\
`nmc`` ``<-`` `[`mode_cluster_state`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md)`(``npost``)`\
[`print`](https://rdrr.io/r/base/print.html)`(``nmc``)`\
`#> <nested mode cluster state: K = 2 modes (Be = 1.98), M = 250 partitions, DL = 13558.3676 nats>`\
`#>  mode   M weight B uncertainty posterior_entropy description_length`\
`#>     1 141  0.564 8   0.1131068          38.15999           6569.625`\
`#>     2 109  0.436 9   0.1579888          48.93175           6809.450`\
\
`# Reconstruct a clean consensus hierarchy from the dominant posterior mode`\
`dom_n`` ``<-`` ``nmc``$``modes``[[``1L``]``]`\
`cons_levels`` ``<-`` `[`mode_consensus`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md)`(``dom_n``)`\
`clean`` ``<-`` `[`list`](https://rdrr.io/r/base/list.html)`(``)`\
`used`` ``<-`` ``NULL`\
`for`` ``(``l`` ``in`` `[`seq_along`](https://rdrr.io/r/base/seq.html)`(``cons_levels``)``)`` ``{`\
`  ``b`` ``<-`` ``cons_levels``[[``l``]``]`\
`  ``if`` ``(``l`` ``>`` ``1L``)`` ``b`` ``<-`` ``b``[``used``]`\
`  ``b``[`[`is.na`](https://rdrr.io/r/base/NA.html)`(``b``)``]`` ``<-`` ``1L`\
`  ``used`` ``<-`` `[`unique`](https://rdrr.io/r/base/unique.html)`(``b``)`\
`  ``clean``[[``l``]``]`` ``<-`` `[`match`](https://rdrr.io/r/base/match.html)`(``b``, ``used``)`\
`  ``if`` ``(`[`length`](https://rdrr.io/r/base/length.html)`(``used``)`` ``==`` ``1L``)`` ``break`\
`}`\
`cons_state`` ``<-`` `[`nested_block_state`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md)`(``g``, ``clean``, degree_corrected ``=`` ``TRUE``)`\
\
`old_par`` ``<-`` ``graphics``::`[`par`](https://rdrr.io/r/graphics/par.html)`(``mfrow ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``1L``, ``2L``)``, mar ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``2``, ``2``, ``3``, ``1``)``)`\
[`plot_hierarchy`](https://diogro.github.io/rgraphtool/reference/plot_hierarchy.md)`(``cons_state``,`\
`               main ``=`` `[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"Dominant mode hierarchy (w = %.2f)"``, ``nmc``$``weights``[``1L``]``)``)`\
[`plot_marginals`](https://diogro.github.io/rgraphtool/reference/plot_marginals.md)`(``g``, ``dom_n``, positions ``=`` ``lesmis``$``positions``, radius ``=`` ``0.045``,`\
`               main ``=`` `[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"Level-0 marginals (unc = %.3f)"``,`\
`                              `[`mode_uncertainty`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md)`(``dom_n``)``)``)`

![](posterior-sampling_files/figure-html/nested-posterior-1.png)

\
`graphics``::`[`par`](https://rdrr.io/r/graphics/par.html)`(``old_par``)`
