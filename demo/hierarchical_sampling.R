# Hierarchical posterior sampling and nested mode characterization
# (graph-tool cookbook, "Sampling from the posterior distribution" - hierarchical
# section). Data are downloaded from Netzschleuder (lesmis).
library(rgraphtool)

lesmis <- netzschleuder_network("lesmis")
g <- lesmis$graph

# 1. Fit nested SBM, equilibrate with nested merge-split MCMC, and sample
set.seed(42)
nstate <- minimize_nested_blockmodel_dl(g)
print(nstate)

npost <- sample_posterior(nstate, nsamples = 400, thin = 5,
                          equilibrate_args = list(wait = 100, nbreaks = 2))

# 2. Posterior distribution of level-0 group count B_0 and hierarchy depth
B0 <- vapply(npost$samples, function(bs) max(bs[[1]]), integer(1))
depth <- vapply(npost$samples, function(bs) {
  top <- match(1L, vapply(bs, max, integer(1)))
  if (is.na(top)) length(bs) else top
}, integer(1))
cat("Nested posterior P(B_0 | A):\n")
print(round(prop.table(table(B0)), 3))
cat("Nested posterior hierarchy depth:\n")
print(round(prop.table(table(depth)), 3))

# 3. Hierarchical posterior modes and consensus hierarchy
nmc <- mode_cluster_state(npost)
print(nmc)

# 4. Visualize the dominant mode's consensus hierarchy and level-0 vertex marginals
dom <- nmc$modes[[1]]
cons_levels <- mode_consensus(dom)
clean <- list()
used <- NULL
for (l in seq_along(cons_levels)) {
  b <- cons_levels[[l]]
  if (l > 1L) b <- b[used]
  b[is.na(b)] <- 1L
  used <- unique(b)
  clean[[l]] <- match(b, used)
  if (length(used) == 1L) break
}
cons_state <- nested_block_state(g, clean, degree_corrected = TRUE)

old_par <- graphics::par(mfrow = c(1, 2), mar = c(2, 2, 3, 1))
on.exit(graphics::par(old_par), add = TRUE)
plot_hierarchy(cons_state, main = sprintf("Dominant mode hierarchy (w = %.2f)",
                                          nmc$weights[1]))
plot_marginals(g, dom, positions = lesmis$positions, radius = 0.045,
               main = sprintf("Level-0 marginals (unc = %.3f)",
                              mode_uncertainty(dom)))
