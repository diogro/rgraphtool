# Weighted stochastic block models and covariate model selection on Les
# Miserables (graph-tool cookbook, "Edge weights and covariates").
# Data are downloaded from Netzschleuder (lesmis: 77 characters, 254 weighted
# co-occurrence edges with integer edge property "value").
library(graphtool)

lesmis <- netzschleuder_network("lesmis")
g <- lesmis$graph

cat("Edge properties available:", paste(graph_edge_properties(g), collapse = ", "), "\n")
w <- graph_edge_property(g, "value")
cat("Edge weight summary (scene co-occurrences):\n")
print(summary(w))

# 1. Compare discrete count covariate models by MAP description length
#    (auto-detected because "value" contains non-negative integers)
set.seed(42)
cmp_disc <- compare_weight_models(
  g,
  weights = "value",
  nested = TRUE,
  method = "map",
  nruns = 3L
)
cat("\nDiscrete count weight models (nested, MAP):\n")
print(cmp_disc)

# 2. Compare continuous and transformed weight models with exact
#    change-of-variables Jacobian adjustment (log-normal vs exponential/normal)
#    using full Bayesian posterior marginal evidence
set.seed(42)
cmp_cont <- compare_weight_models(
  g,
  weights = "value",
  candidates = c("real-exponential", "real-normal", "log-normal"),
  nested = TRUE,
  method = "posterior",
  nsamples = 100L,
  thin = 5L,
  burnin = 50L
)
cat("\nContinuous / transformed weight models (nested, posterior evidence):\n")
print(cmp_cont)

# 3. Posterior sampling and mode clustering on the best discrete weighted model,
#    plus comparison with the unweighted nested SBM
set.seed(42)
best_w <- cmp_disc$best_state
post_w <- sample_posterior(best_w, nsamples = 100L, thin = 5L, burnin = 50L)
modes_w <- mode_cluster_state(post_w$samples, nested = TRUE)
cat("\nPosterior modes for best discrete weighted model (", cmp_disc$best_model, "):\n",
    sep = "")
print(modes_w)

set.seed(42)
unw <- minimize_nested_blockmodel_dl(g)
ov <- partition_overlap(nested_levels(best_w)[[1L]], nested_levels(unw)[[1L]],
                        norm = TRUE)
cat(sprintf("\nLevel-0 partition overlap (weighted %s vs unweighted): %.3f\n",
            cmp_disc$best_model, ov))

# 4. Visualizations: continuous/log-normal evidence comparison and weighted
#    hierarchical community structure
old_par <- graphics::par(mfrow = c(1, 2), mar = c(7, 4, 3, 1))
on.exit(graphics::par(old_par), add = TRUE)

tab <- cmp_cont$summary
graphics::barplot(
  tab$delta,
  names.arg = tab$model,
  las = 2,
  ylab = "Delta(-log P(A, w)) [nats]",
  col = ifelse(tab$model == cmp_cont$best_model, "darkorange", "steelblue"),
  main = "Continuous weight models (0 = best)"
)

graphics::par(mar = c(2, 2, 3, 1))
plot_hierarchy(
  best_w,
  main = sprintf("Les Miserables weighted nested SBM (%s)", cmp_disc$best_model)
)
