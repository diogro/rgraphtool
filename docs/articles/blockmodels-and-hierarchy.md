# Stochastic Block Models & Hierarchical Partitions

This article brings together the `football`, `celegans`, and
`model_selection` cookbook workflows, introducing flat and hierarchical
(nested) stochastic block models (SBMs) and minimum-description-length
(MDL) model selection between degree-corrected (DC) and
non-degree-corrected (NDC) variants.

## 1. Flat SBM inference on the college football network

In a flat SBM, each vertex \\(i\\) is assigned to a group \\(b_i \in
\\{1, \dots, B\\}\\), and we minimize the exact microcanonical
description length \\(\Sigma(A, b)\\) combining the graph likelihood and
the prior entropy of the partition and block edge counts.
[`compare_blockmodels()`](https://diogro.github.io/rgraphtool/reference/minimize_blockmodel_dl.md)
fits both degree-corrected and non-degree-corrected models and computes
the posterior odds ratio between them.

\
`football`` ``<-`` ``load_or_fallback_football``(``)`\
`g`` ``<-`` ``football``$``graph`\
\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``43``)`\
`comparison`` ``<-`` `[`compare_blockmodels`](https://diogro.github.io/rgraphtool/reference/minimize_blockmodel_dl.md)`(``g``)`\
[`print`](https://rdrr.io/r/base/print.html)`(``comparison``$``comparison``)`\
`#>                  model description_length`\
`#> 1 non-degree-corrected           1740.158`\
`#> 2     degree-corrected           1787.537`\
\
`state`` ``<-`` ``comparison``$``best_state`\
[`block_sizes`](https://diogro.github.io/rgraphtool/reference/block_state.md)`(``state``)`\
`#>  1  2  3  4  5  6  7  8  9 10 `\
`#> 14  9 11 12 10 14 12 10 14  9`\
\
`pos`` ``<-`` ``football``$``positions`\
`edges`` ``<-`` `[`graph_edges`](https://diogro.github.io/rgraphtool/reference/graph.md)`(``g``)`\
`groups`` ``<-`` `[`match`](https://rdrr.io/r/base/match.html)`(`[`block_membership`](https://diogro.github.io/rgraphtool/reference/block_state.md)`(``state``)``, `[`sort`](https://rdrr.io/r/base/sort.html)`(`[`unique`](https://rdrr.io/r/base/unique.html)`(`[`block_membership`](https://diogro.github.io/rgraphtool/reference/block_state.md)`(``state``)``)``)``)`\
`colors`` ``<-`` ``grDevices``::`[`hcl.colors`](https://rdrr.io/r/grDevices/palettes.html)`(`[`max`](https://rdrr.io/r/base/Extremes.html)`(``groups``)``, ``"Dynamic"``)`\
\
`old_par`` ``<-`` ``graphics``::`[`par`](https://rdrr.io/r/graphics/par.html)`(``mfrow ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``1L``, ``2L``)``, mar ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``3.5``, ``3.5``, ``2.5``, ``1``)``)`\
`graphics``::`[`plot`](https://rdrr.io/r/graphics/plot.default.html)`(``pos``, type ``=`` ``"n"``, xlab ``=`` ``""``, ylab ``=`` ``""``, axes ``=`` ``FALSE``, asp ``=`` ``1``)`\
`graphics``::`[`segments`](https://rdrr.io/r/graphics/segments.html)`(``pos``[``edges``[``, ``1L``]``, ``1L``]``, ``pos``[``edges``[``, ``1L``]``, ``2L``]``,`\
`                   ``pos``[``edges``[``, ``2L``]``, ``1L``]``, ``pos``[``edges``[``, ``2L``]``, ``2L``]``, col ``=`` ``"grey82"``)`\
`graphics``::`[`points`](https://rdrr.io/r/graphics/points.html)`(``pos``, pch ``=`` ``19L``, col ``=`` ``colors``[``groups``]``, cex ``=`` ``0.85``)`\
`graphics``::`[`title`](https://rdrr.io/r/graphics/title.html)`(``"Football network: inferred blocks"``)`\
\
`counts`` ``<-`` `[`block_edge_counts`](https://diogro.github.io/rgraphtool/reference/block_state.md)`(``state``)`\
`graphics``::`[`image`](https://rdrr.io/r/graphics/image.html)`(`[`seq_len`](https://rdrr.io/r/base/seq.html)`(`[`nrow`](https://rdrr.io/r/base/nrow.html)`(``counts``)``)``, `[`seq_len`](https://rdrr.io/r/base/seq.html)`(`[`ncol`](https://rdrr.io/r/base/nrow.html)`(``counts``)``)``, ``counts``,`\
`                xlab ``=`` ``"Block"``, ylab ``=`` ``"Block"``, main ``=`` ``"Block edge counts"``)`

![](blockmodels-and-hierarchy_files/figure-html/football-fit-1.png)

\
`graphics``::`[`par`](https://rdrr.io/r/graphics/par.html)`(``old_par``)`

## 2. Hierarchical (nested) SBM on *C. elegans*

Flat SBMs have a resolution limit proportional to \\(\sqrt{N}\\) when
\\(B\\) grows. The nested SBM resolves this by modeling the block
multigraph recursively as another SBM at higher levels \\(\ell = 1, 2,
\dots, L\\).

\
`celegans`` ``<-`` ``load_or_fallback_celegans``(``)`\
`g_ce`` ``<-`` ``celegans``$``graph`\
\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``1``)`\
`nstate`` ``<-`` `[`minimize_nested_blockmodel_dl`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md)`(``g_ce``)`\
[`print`](https://rdrr.io/r/base/print.html)`(``nstate``)`\
`#> <nested block state: degree-corrected, 10 levels, description length 8596.5625 nats>`\
`#> l: 0, N: 297, B: 22`\
`#> l: 1, N: 22, B: 4`\
`#> l: 2, N: 4, B: 1`\
`#> l: 3, N: 1, B: 1`\
[`nested_summary`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md)`(``nstate``)`\
`#>   level   N  B description_length`\
`#> 1     0 297 22         7538.23546`\
`#> 2     1  22  4          968.29937`\
`#> 3     2   4  1           90.02768`\
`#> 4     3   1  1            0.00000`\
\
`# Project a single vertex's group membership across every level of the hierarchy`\
`v_id`` ``<-`` `[`min`](https://rdrr.io/r/base/Extremes.html)`(``47L``, `[`graph_num_vertices`](https://diogro.github.io/rgraphtool/reference/graph.md)`(``g_ce``)``)`\
[`sapply`](https://rdrr.io/r/base/lapply.html)`(`[`seq_len`](https://rdrr.io/r/base/seq.html)`(`[`nrow`](https://rdrr.io/r/base/nrow.html)`(`[`nested_summary`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md)`(``nstate``)``)``)`` ``-`` ``1L``,`\
`       ``function``(``l``)`` `[`nested_project`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md)`(``nstate``, ``l``)``[``v_id``]``)`\
`#> [1] 9 1 1 1`\
\
`# Inspect level-0 group sizes as a flat block state`\
`level0`` ``<-`` `[`nested_level_state`](https://diogro.github.io/rgraphtool/reference/minimize_nested_blockmodel_dl.md)`(``nstate``, ``0L``)`\
[`block_sizes`](https://diogro.github.io/rgraphtool/reference/block_state.md)`(``level0``)`\
`#>  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 `\
`#> 12 10  7  5  9 11 21 20 24 20 12  7  3 15 14 14 15 15 23 19  9 12`\
\
[`plot_hierarchy`](https://diogro.github.io/rgraphtool/reference/plot_hierarchy.md)`(``nstate``, main ``=`` ``"C. elegans: nested SBM"``)`

![](blockmodels-and-hierarchy_files/figure-html/celegans-nested-1.png)

## 3. Degree-corrected vs. non-degree-corrected nested model selection

Using `compare_blockmodels(..., nested = TRUE)` evaluates whether
empirical degree heterogeneity warrants the degree-corrected
microcanonical prior. On networks with broad degree distributions (such
as political blogs), the DC nested SBM achieves a decisively lower
description length; on nearly regular networks (such as college
football), the simpler NDC model is often preferred.

\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``1``)`\
`cmp_fb`` ``<-`` `[`compare_blockmodels`](https://diogro.github.io/rgraphtool/reference/minimize_blockmodel_dl.md)`(``football``$``graph``, nested ``=`` ``TRUE``, nruns ``=`` ``3L``)`\
[`print`](https://rdrr.io/r/base/print.html)`(``cmp_fb``)`\
`#> <nested block model comparison>`\
`#> non-degree-corrected:  description length    1733.5257 nats`\
`#> degree-corrected:      description length    1782.0200 nats`\
`#> ln Lambda (NDC vs DC): 48.4943, about 10^21 in favour of the non-degree-corrected model`\
`cmp_fb``$``comparison`\
`#>                  model description_length`\
`#> 1 non-degree-corrected           1733.526`\
`#> 2     degree-corrected           1782.020`

When network access to Netzschleuder is available, the same comparison
runs on `polblogs` (`1,490` blogs, `19,090` directed links):

\
`polblogs`` ``<-`` `[`netzschleuder_network`](https://diogro.github.io/rgraphtool/reference/netzschleuder_network.md)`(``"polblogs"``)`\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``1``)`\
[`compare_blockmodels`](https://diogro.github.io/rgraphtool/reference/minimize_blockmodel_dl.md)`(``polblogs``$``graph``, nested ``=`` ``TRUE``)`
