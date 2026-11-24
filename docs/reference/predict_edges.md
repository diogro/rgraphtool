# Bayesian missing and spurious edge prediction and evaluation

Computes model-averaged Bayesian posterior edge probabilities for
candidate missing or spurious links using flat/nested SBMs or
uncertain-network models, and evaluates binary link-prediction
performance with exact tie-aware ROC-AUC, PR-AUC, Brier score, and
ROC/PR curves.

## Usage

``` r
predict_edges(x, candidates, nsamples = 100L, burnin = 50L, thin = 1L,
  nested = TRUE, degree_corrected = TRUE,
  degree_prior = c("distributed", "uniform"),
  self_loops = FALSE, epsilon = 1e-8)

evaluate_edge_prediction(scores, labels, threshold = 0.5)
```

## Arguments

- x:

  A `rgraphtool_graph`, a fitted `rgraphtool_block_state` or
  `rgraphtool_nested_block_state`, or a
  `rgraphtool_uncertain_block_state`.

- candidates:

  Two-column integer matrix or data.frame of 1-based candidate vertex
  pairs `(u, v)` to score.

- nsamples, burnin, thin:

  Posterior MCMC sample count, burn-in sweeps, and thinning interval.

- nested, degree_corrected, degree_prior:

  SBM prior options used when `x` is a raw graph.

- self_loops:

  Whether self-loops are valid candidate dyads.

- epsilon:

  Multiplicity summation convergence tolerance.

- scores:

  Numeric vector of predicted edge scores/probabilities, or a
  `rgraphtool_edge_prediction` object.

- labels:

  Binary ground-truth labels (`0`/`1` or logical) with at least one
  positive and one negative example.

- threshold:

  Classification threshold for precision, recall, F1, and accuracy.

## Value

`predict_edges()` returns a `rgraphtool_edge_prediction` containing a
`predictions` data.frame (`from`, `to`, `prob`, `log_prob`, and
`observed` or `conditional_prob`).

`evaluate_edge_prediction()` returns a `rgraphtool_edge_evaluation`
containing `roc_auc`, `pr_auc`, `precision`, `recall`, `f1`, `accuracy`,
`brier`, `confusion`, and full `roc` and `pr` curve data.frames.

## References

Peixoto, T. P. (2018). Reconstructing networks with unknown and
heterogeneous errors. *Physical Review X* 8, 041011.
doi:10.1103/PhysRevX.8.041011

## See also

[`measured_block_state`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md),
[`edge_probabilities`](https://diogro.github.io/rgraphtool/reference/uncertain_block_state.md),
`demo("uncertain_reconstruction")`

## Examples

``` r
g <- graph_create(12, directed = FALSE)
e <- rbind(t(combn(1:6, 2)), t(combn(7:12, 2)))
graph_add_edges(g, e[-1, , drop = FALSE])
cands <- rbind(e[1, ], c(1L, 12L))
set.seed(1)
pred <- predict_edges(g, cands, nsamples = 20L, burnin = 10L)
evaluate_edge_prediction(pred, c(1L, 0L))
#> <rgraphtool_edge_evaluation>
#>   ROC-AUC: 1.0000 | PR-AUC: 1.0000 | Brier: 0.0042
#>   At threshold 0.50 -> Precision: 1.0000 | Recall: 1.0000 | F1: 1.0000 | Acc: 1.0000
```
