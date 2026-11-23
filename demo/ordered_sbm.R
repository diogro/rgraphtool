# Ordered community detection with the flat and nested ranked SBM on the
# C. elegans neural network (graph-tool cookbook, "Ordered community structure";
# Peixoto 2022).
# Data are downloaded from Netzschleuder (celegansneural: 297 neurons,
# 2,359 directed synaptic connections).
library(graphtool)

celegans <- netzschleuder_network("celegansneural")
g <- celegans$graph

cat(sprintf("C. elegans neural network: N = %d neurons, E = %d directed edges\n",
            graph_num_vertices(g), graph_num_edges(g)))

# 1. Compare unordered vs. ordered (ranked) flat and nested SBMs by MAP
#    minimum description length
set.seed(42)
cmp_map <- compare_ordered_models(
  g,
  models = c("unordered_flat", "unordered_nested",
             "ranked_flat", "ranked_nested"),
  method = "map",
  nruns = 2L
)
cat("\nUnordered vs. ordered SBM comparison (MAP):\n")
print(cmp_map)

# 2. Full Bayesian posterior marginal evidence comparing unordered vs. ranked
#    nested SBMs
set.seed(42)
cmp_post <- compare_ordered_models(
  g,
  models = c("unordered_nested", "ranked_nested"),
  method = "posterior",
  nsamples = 60L,
  thin = 5L,
  burnin = 40L
)
cat("\nUnordered vs. ordered nested SBM (posterior marginal evidence):\n")
print(cmp_post)

# 3. Inspect the inferred feed-forward ordering of the nested ranked SBM
rk_nest <- cmp_map$states$ranked_nested
ec <- ranked_edge_counts(rk_nest)
ec_tot <- ec$upstream + ec$lateral + ec$downstream
cat(sprintf(
  "\nNested ranked SBM edge flow: %d upstream (%.1f%%), %d lateral (%.1f%%), %d downstream (%.1f%%)\n",
  ec$upstream, 100 * ec$upstream / ec_tot,
  ec$lateral, 100 * ec$lateral / ec_tot,
  ec$downstream, 100 * ec$downstream / ec_tot
))
cat("\nHierarchy summary:\n")
print(nested_summary(rk_nest))

# 4. Visualizations: MAP description length comparison and rank-ordered layout
old_par <- graphics::par(mfrow = c(1, 2), mar = c(8, 4, 3, 1))
on.exit(graphics::par(old_par), add = TRUE)

tab <- cmp_map$summary
graphics::barplot(
  tab$delta,
  names.arg = tab$model,
  las = 2,
  ylab = "Delta(DL) [nats]",
  col = ifelse(tab$model == cmp_map$best_model, "darkorange", "steelblue"),
  main = "Unordered vs. ordered SBM (0 = best)"
)

graphics::par(mar = c(1, 1, 3, 1))
plot_ranked_hierarchy(
  rk_nest,
  vertex_cex = 0.65,
  main = "C. elegans ordered community hierarchy"
)
