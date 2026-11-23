# Layered stochastic block models and edge-covariate vs. independent-layers
# model selection on Les Miserables (graph-tool cookbook, "Layered networks").
# Data are downloaded from Netzschleuder (lesmis: 77 characters, 254 edges
# stratified across discrete scene co-occurrence layers "value").
library(graphtool)

lesmis <- netzschleuder_network("lesmis")
g <- lesmis$graph

w <- graph_edge_property(g, "value")
cat(sprintf("Les Miserables: N = %d vertices, E = %d edges, C = %d discrete layers\n",
            graph_num_vertices(g), graph_num_edges(g), length(unique(w))))
print(table(layer = w))

# 1. Compare edge-covariate vs. independent-layers models (DC and NDC) by MAP
#    minimum description length (3 restarts each)
set.seed(42)
cmp_map <- compare_layered_models(
  g,
  layers = "value",
  models = c("edge_covariates_dc", "edge_covariates_ndc",
             "independent_layers_dc", "independent_layers_ndc"),
  nested = TRUE,
  method = "map",
  nruns = 3L
)
cat("\nLayered model comparison (nested, MAP):\n")
print(cmp_map)

# 2. Full Bayesian posterior marginal evidence (-log P(A, l) = <Sigma> - H(b))
#    comparing edge-covariates vs. independent-layers
set.seed(42)
cmp_post <- compare_layered_models(
  g,
  layers = "value",
  models = c("edge_covariates", "independent_layers"),
  nested = TRUE,
  method = "posterior",
  nsamples = 100L,
  thin = 5L,
  burnin = 50L
)
cat("\nLayered model comparison (nested, posterior marginal evidence):\n")
print(cmp_post)

# 3. Posterior mode clustering on the winning layered model and comparison
#    with the unweighted nested SBM
set.seed(42)
best_l <- cmp_post$best_state
post_l <- sample_posterior(best_l, nsamples = 100L, thin = 5L, burnin = 50L)
modes_l <- mode_cluster_state(post_l$samples, nested = TRUE)
cat("\nPosterior modes for best layered model (", cmp_post$best_model, "):\n",
    sep = "")
print(modes_l)

set.seed(42)
unw <- minimize_nested_blockmodel_dl(g)
ov <- partition_overlap(nested_levels(best_l)[[1L]], nested_levels(unw)[[1L]],
                        norm = TRUE)
cat(sprintf("\nLevel-0 partition overlap (%s vs unweighted): %.3f\n",
            cmp_post$best_model, ov))

# 4. Visualizations: MAP description length across layered model variants and
#    hierarchical community structure of the winning layered model
old_par <- graphics::par(mfrow = c(1, 2), mar = c(9, 4, 3, 1))
on.exit(graphics::par(old_par), add = TRUE)

tab <- cmp_map$summary
graphics::barplot(
  tab$delta,
  names.arg = tab$model,
  las = 2,
  ylab = "Delta(DL) [nats]",
  col = ifelse(tab$model == cmp_map$best_model, "darkorange", "steelblue"),
  main = "Layered model variants (0 = best)"
)

graphics::par(mar = c(2, 2, 3, 1))
plot_hierarchy(
  best_l,
  main = sprintf("Les Miserables layered SBM (%s)", cmp_post$best_model)
)
