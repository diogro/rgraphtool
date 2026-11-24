# Fungus-Tree Ecological Network Case Study (Vacher et al. 2008; GrossSBM/sbm):
# Bipartite SBM on the 154 x 51 fungus-tree interaction matrix, simple &
# covariate-weighted SBMs on the 51 x 51 tree-tree network, and posterior MCMC.

library(rgraphtool)

set.seed(42)

# 1. Load the bundled Vacher et al. (2008) fungus-tree dataset
ft <- fungus_tree_network()
cat(sprintf("Loaded fungus-tree network: %d fungi x %d trees (%d interactions)\n",
            nrow(ft$fungus_tree), ncol(ft$fungus_tree), sum(ft$fungus_tree)))
print(ft$graph, max_edges = 8L)

# 2. Fit microcanonical Bipartite SBM directly to the 154 x 51 incidence matrix
fit_bip <- minimize_bipartite_blockmodel_dl(
  ft$fungus_tree,
  dim_labels = c("fungus", "tree"),
  nruns = 3L
)
print(fit_bip)
cat("\nBlock interaction densities (fungus blocks x tree blocks):\n")
print(round(bipartite_block_matrix(fit_bip, normalize = "density"), 3))

# 3. Compare flat vs. nested and degree-corrected vs. non-degree-corrected models
cmp_bip <- compare_bipartite_models(
  ft$fungus_tree,
  dim_labels = c("fungus", "tree"),
  nruns = 2L
)
print(cmp_bip)

# 4. Posterior MCMC sampling over bipartite partitions
post_bip <- sample_posterior(fit_bip, nsamples = 40L, burnin = 20L, thin = 2L)
cat(sprintf("Posterior mean DL: %.2f nats (sd %.2f)\n\n",
            mean(post_bip$entropies), stats::sd(post_bip$entropies)))

# 5. Unipartite tree-tree network: binary vs. Poisson count + covariates
tree_mat <- ft$tree_tree
diag(tree_mat) <- 0
idx_ut <- which(upper.tri(tree_mat) & tree_mat > 0, arr.ind = TRUE)
g_tree <- graph_create(nrow(tree_mat), directed = FALSE)
graph_set_vertex_property(g_tree, "name", ft$tree_names)
graph_add_edges(
  g_tree,
  idx_ut,
  properties = list(
    shared_fungi = tree_mat[idx_ut],
    genetic_dist = ft$covar_tree$genetic_dist[idx_ut],
    taxonomic_dist = ft$covar_tree$taxonomic_dist[idx_ut],
    geographic_dist = ft$covar_tree$geographic_dist[idx_ut]
  )
)

fit_tree_poi <- minimize_weighted_blockmodel_dl(
  g_tree, weights = "shared_fungi", types = "discrete-poisson"
)
fit_tree_cov <- minimize_weighted_blockmodel_dl(
  g_tree,
  weights = c("shared_fungi", "taxonomic_dist"),
  types = c("discrete-poisson", "real-exponential")
)
cat("Tree-tree weighted SBM (Poisson shared fungi):\n")
print(fit_tree_poi)
cat("Tree-tree multi-covariate SBM (shared fungi + taxonomic distance):\n")
print(fit_tree_cov)

# 6. Visualize block-reordered fungus-tree matrix and mesoscale block network
old_par <- graphics::par(mfrow = c(1L, 2L), mar = c(4, 4, 3, 1))
on.exit(graphics::par(old_par), add = TRUE)
plot_bipartite(fit_bip, type = "matrix", main = "Fungus-Tree Reordered Matrix")
plot_bipartite(fit_bip, type = "meso", main = "Fungus-Tree Mesoscale Blocks")
