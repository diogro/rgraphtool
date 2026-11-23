# Posterior sampling, vertex marginals, and posterior modes on Les Miserables
# (graph-tool cookbook, "Sampling from the posterior distribution").
# Data are downloaded from Netzschleuder (lesmis: 77 characters, 254 edges).
library(rgraphtool)

lesmis <- netzschleuder_network("lesmis")
g <- lesmis$graph

# 1. Start from a multilevel point estimate, then equilibrate merge-split MCMC
set.seed(42)
state <- minimize_blockmodel_dl(g)
eq <- mcmc_equilibrate(state, wait = 100, nbreaks = 2,
                       mcmc_args = list(niter = 10))
cat(sprintf("Equilibrated flat DC SBM in %d steps (DL = %.2f nats)\n",
            eq$niter, eq$entropy))

# 2. Collect posterior samples (every 10 sweeps)
post <- sample_posterior(eq$state, nsamples = 500, thin = 10, burnin = 0)

# Posterior distribution of the number of groups B
B_samples <- vapply(post$samples, function(b) max(b), integer(1))
B_dist <- prop.table(table(B_samples))
cat("Posterior P(B | A):\n")
print(round(B_dist, 3))

# 3. Single-mode alignment, vertex marginals, and consensus partition
pm <- partition_mode_state(post)
print(pm)

# 4. Posterior mode clustering (detecting multimodal structure)
mc <- mode_cluster_state(post)
print(mc)

# 5. Visualizations: P(B | A) and pie-chart vertex marginals of the dominant mode
old_par <- graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
on.exit(graphics::par(old_par), add = TRUE)

graphics::barplot(B_dist, xlab = "Number of groups B", ylab = "Posterior P(B | A)",
                  col = "steelblue", main = "Les Miserables: P(B | A)")

dom <- mc$modes[[1]]
plot_marginals(g, dom, positions = lesmis$positions, radius = 0.045,
               main = sprintf("Dominant mode (w = %.2f, unc = %.3f)",
                              mc$weights[1], mode_uncertainty(dom)))
