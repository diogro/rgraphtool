# Assortative community detection with the planted-partition SBM and model
# comparison against the general SBM on Les Miserables (graph-tool cookbook,
# "Assortative community structure"; Zhang & Peixoto 2020).
# Data are downloaded from Netzschleuder (lesmis: 77 characters, 254 edges).
library(rgraphtool)

lesmis <- netzschleuder_network("lesmis")
g <- lesmis$graph

cat(sprintf("Les Miserables: N = %d vertices, E = %d edges\n",
            graph_num_vertices(g), graph_num_edges(g)))

# 1. Fit non-uniform and uniform planted-partition models
set.seed(42)
pp_nu <- minimize_pp_blockmodel_dl(g, uniform = FALSE)
pp_u  <- minimize_pp_blockmodel_dl(g, uniform = TRUE)

cat("\nNon-uniform planted-partition fit:\n")
print(pp_nu)
print(pp_assortativity(pp_nu))

cat("\nUniform planted-partition fit:\n")
print(pp_u)
print(pp_assortativity(pp_u))

# 2. Compare assortative planted-partition vs. general flat and nested SBMs
#    by MAP minimum description length
set.seed(42)
cmp_map <- compare_assortative_models(
  g,
  models = c("pp", "pp_uniform", "sbm_flat", "sbm_nested"),
  method = "map",
  nruns = 3L
)
cat("\nAssortative vs. general SBM comparison (MAP):\n")
print(cmp_map)

# 3. Full Bayesian posterior marginal evidence (-log P(A) = <Sigma> - H(b))
set.seed(42)
cmp_post <- compare_assortative_models(
  g,
  models = c("pp", "pp_uniform", "sbm_flat", "sbm_nested"),
  method = "posterior",
  nsamples = 100L,
  thin = 5L,
  burnin = 50L
)
cat("\nAssortative vs. general SBM comparison (posterior marginal evidence):\n")
print(cmp_post)

# 4. Posterior sampling and consensus mode of the planted-partition model
set.seed(42)
post_pp <- sample_posterior(pp_nu, nsamples = 100L, thin = 5L, burnin = 50L)
mode_pp <- partition_mode_state(post_pp$samples)
cat("\nPlanted-partition posterior consensus mode:\n")
print(mode_pp)

# 5. Visualizations: posterior model comparison and vertex membership marginals
old_par <- graphics::par(mfrow = c(1, 2), mar = c(7, 4, 3, 1))
on.exit(graphics::par(old_par), add = TRUE)

tab <- cmp_post$summary
graphics::barplot(
  tab$delta,
  names.arg = tab$model,
  las = 2,
  ylab = "Delta(-log P(A)) [nats]",
  col = ifelse(tab$model == cmp_post$best_model, "darkorange", "steelblue"),
  main = "Assortative vs. general SBM (0 = best)"
)

graphics::par(mar = c(2, 2, 3, 1))
plot_marginals(
  g,
  mode_pp,
  main = sprintf("Les Miserables planted-partition marginals (B = %d)", pp_nu$B)
)
