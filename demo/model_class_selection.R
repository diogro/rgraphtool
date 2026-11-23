# Model class selection across flat, nested, degree-corrected, and overlapping
# stochastic block models on Les Miserables (graph-tool cookbook, "Model class
# selection"). Data are downloaded from Netzschleuder (lesmis: 77 characters,
# 254 edges).
library(graphtool)

lesmis <- netzschleuder_network("lesmis")
g <- lesmis$graph

# 1. Compare all 8 model classes by MAP minimum description length (3 restarts each)
set.seed(42)
cmp_map <- compare_model_classes(g, method = "map", nruns = 3L)
cat("MAP minimum description length comparison (all 8 model classes):\n")
print(cmp_map)

# 2. Full Bayesian posterior model evidence (-log P(A) = <Sigma> - H(b))
#    comparing flat vs nested, DC vs NDC, and overlapping vs non-overlapping
set.seed(42)
cmp_post <- compare_model_classes(
  g,
  models = c("flat_ndc", "flat_dc", "nested_ndc", "nested_dc",
             "overlap_flat_dc", "overlap_nested_dc"),
  method = "posterior",
  nsamples = 200L,
  thin = 5L,
  burnin = 50L
)
cat("\nPosterior marginal evidence comparison:\n")
print(cmp_post)

cat("\nPairwise log Bayes factor matrix (positive favors row over column):\n")
print(round(cmp_post$log_bayes_factors, 2))

# 3. Visualizations: relative -log P(A) across model classes and overlapping
#    mixed memberships on Les Miserables
old_par <- graphics::par(mfrow = c(1, 2), mar = c(8, 4, 3, 1))
on.exit(graphics::par(old_par), add = TRUE)

tab <- cmp_post$summary
graphics::barplot(
  tab$delta,
  names.arg = tab$model,
  las = 2,
  ylab = "Delta(-log P(A)) [nats]",
  col = ifelse(tab$model == cmp_post$best_model, "darkorange", "steelblue"),
  main = "Model class evidence (0 = best)"
)

ov_st <- cmp_post$states$overlap_flat_dc
ov_mat <- overlap_groups(ov_st)$matrix
n_ov <- sum(rowSums(ov_mat > 0L) > 1L)
graphics::par(mar = c(2, 2, 3, 1))
plot_marginals(
  g,
  ov_mat,
  positions = lesmis$positions,
  radius = 0.045,
  main = sprintf("Overlap flat DC (B = %d, %d overlapping)",
                 overlap_groups(ov_st)$B, n_ov)
)
