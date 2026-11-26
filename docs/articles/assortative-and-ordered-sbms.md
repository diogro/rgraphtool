# Assortative & Ordered (Ranked) Communities

This article combines the `assortative_sbm` and `ordered_sbm` demos,
covering two specialized microcanonical SBM priors: the
**planted-partition SBM** for assortative community detection (Zhang &
Peixoto 2020) and the **ranked (ordered) SBM** for discovering directed
feed-forward hierarchies (Peixoto 2022).

## 1. Assortative community structure (planted-partition SBM)

The planted-partition model restricts block interactions to within-group
vs. between-group densities, yielding a sharper prior whenever community
structure is predominantly assortative.
[`compare_assortative_models()`](https://diogro.github.io/rgraphtool/reference/pp_block_state.md)
tests non-uniform and uniform planted-partition models directly against
general flat and nested SBMs.

\
`lesmis`` ``<-`` ``load_or_fallback_lesmis``(``)`\
`g_lm`` ``<-`` ``lesmis``$``graph`\
\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``42``)`\
`pp_nu`` ``<-`` `[`minimize_pp_blockmodel_dl`](https://diogro.github.io/rgraphtool/reference/pp_block_state.md)`(``g_lm``, uniform ``=`` ``FALSE``)`\
`pp_u``  ``<-`` `[`minimize_pp_blockmodel_dl`](https://diogro.github.io/rgraphtool/reference/pp_block_state.md)`(``g_lm``, uniform ``=`` ``TRUE``)`\
[`print`](https://rdrr.io/r/base/print.html)`(``pp_nu``)`\
`#> <planted-partition block state: degree-corrected, non-uniform assortativity, B = 5 (Be = 4.60), E_in/E = 76.8%, DL = 719.4225 nats>`\
[`pp_assortativity`](https://diogro.github.io/rgraphtool/reference/pp_block_state.md)`(``pp_nu``)`\
`#> $E_in`\
`#> [1] 195`\
`#> `\
`#> $E_out`\
`#> [1] 59`\
`#> `\
`#> $E_total`\
`#> [1] 254`\
`#> `\
`#> $fraction_in`\
`#> [1] 0.7677165`\
`#> `\
`#> $B`\
`#> [1] 5`\
`#> `\
`#> $Be`\
`#> [1] 4.598904`\
`#> `\
`#> $err`\
`#> [1] 43 32 29 16 75`\
\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``42``)`\
`cmp_pp`` ``<-`` `[`compare_assortative_models`](https://diogro.github.io/rgraphtool/reference/pp_block_state.md)`(`\
`  ``g_lm``,`\
`  models ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"pp"``, ``"pp_uniform"``, ``"sbm_flat"``, ``"sbm_nested"``)``,`\
`  method ``=`` ``"posterior"``,`\
`  nsamples ``=`` ``80L``,`\
`  thin ``=`` ``5L``,`\
`  burnin ``=`` ``40L`\
`)`\
[`print`](https://rdrr.io/r/base/print.html)`(``cmp_pp``)`\
`#> <assortative vs. general SBM comparison (DC, posterior/variational)>`\
`#>       model B description_length posterior_entropy neg_log_evidence  delta`\
`#>  sbm_nested 6            702.924            29.843          688.607  0.000`\
`#>    sbm_flat 6            700.726            18.937          691.990  3.383`\
`#>          pp 5            719.422            36.175          702.923 14.316`\
`#>  pp_uniform 4            725.074            16.638          712.270 23.664`\
`#>  posterior_weight`\
`#>          9.67e-01`\
`#>          3.28e-02`\
`#>          5.86e-07`\
`#>          5.11e-11`\
`#> Best model: sbm_nested`\
\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``42``)`\
`post_pp`` ``<-`` `[`sample_posterior`](https://diogro.github.io/rgraphtool/reference/mcmc_sweep.md)`(``pp_nu``, nsamples ``=`` ``80L``, thin ``=`` ``5L``, burnin ``=`` ``40L``)`\
`mode_pp`` ``<-`` `[`partition_mode_state`](https://diogro.github.io/rgraphtool/reference/partition_mode_state.md)`(``post_pp``$``samples``)`\
\
`old_par`` ``<-`` ``graphics``::`[`par`](https://rdrr.io/r/graphics/par.html)`(``mfrow ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``1L``, ``2L``)``, mar ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``7``, ``4``, ``3``, ``1``)``)`\
`tab`` ``<-`` ``cmp_pp``$``summary`\
`graphics``::`[`barplot`](https://rdrr.io/r/graphics/barplot.html)`(`\
`  ``tab``$``delta``,`\
`  names.arg ``=`` ``tab``$``model``,`\
`  las ``=`` ``2``,`\
`  ylab ``=`` ``"Delta(-log P(A)) [nats]"``,`\
`  col ``=`` `[`ifelse`](https://rdrr.io/r/base/ifelse.html)`(``tab``$``model`` ``==`` ``cmp_pp``$``best_model``, ``"darkorange"``, ``"steelblue"``)``,`\
`  main ``=`` ``"Assortative vs. general SBM (0 = best)"`\
`)`\
`graphics``::`[`par`](https://rdrr.io/r/graphics/par.html)`(``mar ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``2``, ``2``, ``3``, ``1``)``)`\
[`plot_marginals`](https://diogro.github.io/rgraphtool/reference/plot_marginals.md)`(`\
`  ``g_lm``, ``mode_pp``,`\
`  main ``=`` `[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(``"Planted-partition marginals (B = %d)"``, ``pp_nu``$``B``)`\
`)`

![](assortative-and-ordered-sbms_files/figure-html/assortative-sbm-1.png)

\
`graphics``::`[`par`](https://rdrr.io/r/graphics/par.html)`(``old_par``)`

## 2. Ordered (ranked) communities in directed networks

For directed networks with directional flow (such as neural circuits,
food webs, or citation hierarchies),
[`ranked_block_state()`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md)
and
[`nested_ranked_block_state()`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md)
jointly infer group memberships and a latent rank ordering \\(u_r\\)
separating upstream (\\(u_i \< u_j\\)), lateral (\\(u_i = u_j\\)), and
downstream (\\(u_i \> u_j\\)) edges.

\
`celegans`` ``<-`` ``load_or_fallback_celegans``(``)`\
`g_ce`` ``<-`` ``celegans``$``graph`\
\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``42``)`\
`cmp_ord`` ``<-`` `[`compare_ordered_models`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md)`(`\
`  ``g_ce``,`\
`  models ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"unordered_flat"``, ``"unordered_nested"``,`\
`             ``"ranked_flat"``, ``"ranked_nested"``)``,`\
`  method ``=`` ``"map"``,`\
`  nruns ``=`` ``2L`\
`)`\
[`print`](https://rdrr.io/r/base/print.html)`(``cmp_ord``)`\
`#> <ordered (ranked) vs. unordered SBM comparison (DC, MAP)>`\
`#>             model  B description_length   delta posterior_weight`\
`#>  unordered_nested 25           8540.092   0.000         1.00e+00`\
`#>     ranked_nested 20           8548.436   8.344         2.38e-04`\
`#>       ranked_flat 13           8661.216 121.125         2.49e-53`\
`#>    unordered_flat  9           8718.870 178.779         2.28e-78`\
`#> Best model: unordered_nested`\
\
`rk_nest`` ``<-`` ``cmp_ord``$``states``$``ranked_nested`\
`ec`` ``<-`` `[`ranked_edge_counts`](https://diogro.github.io/rgraphtool/reference/ranked_block_state.md)`(``rk_nest``)`\
`ec_tot`` ``<-`` ``ec``$``upstream`` ``+`` ``ec``$``lateral`` ``+`` ``ec``$``downstream`\
[`cat`](https://rdrr.io/r/base/cat.html)`(`[`sprintf`](https://rdrr.io/r/base/sprintf.html)`(`\
`  ``"Nested ranked SBM flow: %d upstream (%.1f%%), %d lateral (%.1f%%), %d downstream (%.1f%%)\n"``,`\
`  ``ec``$``upstream``, ``100`` ``*`` ``ec``$``upstream`` ``/`` ``ec_tot``,`\
`  ``ec``$``lateral``, ``100`` ``*`` ``ec``$``lateral`` ``/`` ``ec_tot``,`\
`  ``ec``$``downstream``, ``100`` ``*`` ``ec``$``downstream`` ``/`` ``ec_tot`\
`)``)`\
`#> Nested ranked SBM flow: 1153 upstream (48.9%), 1125 lateral (47.7%), 81 downstream (3.4%)`\
\
`old_par`` ``<-`` ``graphics``::`[`par`](https://rdrr.io/r/graphics/par.html)`(``mfrow ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``1L``, ``2L``)``, mar ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``8``, ``4``, ``3``, ``1``)``)`\
`tab_o`` ``<-`` ``cmp_ord``$``summary`\
`graphics``::`[`barplot`](https://rdrr.io/r/graphics/barplot.html)`(`\
`  ``tab_o``$``delta``,`\
`  names.arg ``=`` ``tab_o``$``model``,`\
`  las ``=`` ``2``,`\
`  ylab ``=`` ``"Delta(DL) [nats]"``,`\
`  col ``=`` `[`ifelse`](https://rdrr.io/r/base/ifelse.html)`(``tab_o``$``model`` ``==`` ``cmp_ord``$``best_model``, ``"darkorange"``, ``"steelblue"``)``,`\
`  main ``=`` ``"Unordered vs. ordered SBM (0 = best)"`\
`)`\
`graphics``::`[`par`](https://rdrr.io/r/graphics/par.html)`(``mar ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``1``, ``1``, ``3``, ``1``)``)`\
[`plot_ranked_hierarchy`](https://diogro.github.io/rgraphtool/reference/plot_hierarchy.md)`(`\
`  ``rk_nest``,`\
`  vertex_cex ``=`` ``0.65``,`\
`  main ``=`` ``"C. elegans ordered community hierarchy"`\
`)`

![](assortative-and-ordered-sbms_files/figure-html/ordered-sbm-1.png)

\
`graphics``::`[`par`](https://rdrr.io/r/graphics/par.html)`(``old_par``)`
