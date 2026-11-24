# Bipartite Stochastic Block Models (biSBM): incidence matrices, weighted
# bipartite networks, hierarchical inference, MCMC posterior sampling, and
# matrix/mesoscale visualizations.

library(rgraphtool)

set.seed(42)

# 1. Synthetic plant-pollinator incidence matrix (3 plant x 2 pollinator blocks)
n_plants <- 30L
n_pollinators <- 24L
true_p <- rep(1:3, each = 10L)
true_a <- rep(1:2, each = 12L)

prob_mat <- rbind(
  c(0.75, 0.05),
  c(0.05, 0.75),
  c(0.45, 0.45)
)

M <- matrix(0L, nrow = n_plants, ncol = n_pollinators,
            dimnames = list(plant = paste0("P", seq_len(n_plants)),
                            pollinator = paste0("A", seq_len(n_pollinators))))
for (i in seq_len(n_plants)) {
  for (j in seq_len(n_pollinators)) {
    M[i, j] <- stats::rbinom(1L, 1L, prob_mat[true_p[i], true_a[j]])
  }
}

# Fit directly from the incidence matrix
fit <- minimize_bipartite_blockmodel_dl(M, nruns = 3L)
print(fit)
print(round(bipartite_block_matrix(fit, normalize = "density"), 3))

# 2. Model comparison (DC vs NDC, flat vs nested)
cmp <- compare_bipartite_models(M, nruns = 2L)
print(cmp)

# 3. Posterior MCMC sampling
post <- sample_posterior(fit, nsamples = 40L, burnin = 20L, thin = 2L)
cat(sprintf("Posterior mean DL: %.2f (sd %.2f)\n",
            mean(post$entropies), stats::sd(post$entropies)))

# 4. Visualize reordered incidence matrix and mesoscale block diagram
old_par <- graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
on.exit(graphics::par(old_par), add = TRUE)
plot_bipartite(fit, type = "matrix")
plot_bipartite(fit, type = "meso")
