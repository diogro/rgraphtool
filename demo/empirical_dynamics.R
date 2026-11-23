# Cookbook demo: Empirical network reconstruction from financial time series
# (multivariate Gaussian partial-correlation network via PseudoNormalBlockState).
#
# Preprocessing requirements for PseudoNormalBlockState on financial prices:
#   1. Convert price levels P_i(t) to stationary log-returns:
#      r_i(t) = log(P_i(t+1) / P_i(t)).
#   2. Standardize each series to zero mean and unit variance so inferred
#      precision off-diagonals and L1/MDL penalties are scale-invariant across
#      indices or stocks (fix_mean = TRUE also enforces centering).
#   3. Use positive = TRUE so the reconstructed precision matrix Omega is
#      strictly diagonally dominant and positive-definite.

library(rgraphtool)
set.seed(42)

cat("=== 1. Empirical European stock indices (datasets::EuStockMarkets) ===\n")
prices <- as.matrix(datasets::EuStockMarkets)
# Step 1: log-returns (T-1 x 4) -> transpose to N x M (4 indices x 1859 days)
log_ret <- diff(log(prices))
# Step 2: standardize each index to zero mean and unit variance
z_ret <- t(scale(log_ret))
idx_names <- colnames(prices)
cat(sprintf("Input matrix: %d assets (%s) x %d trading days\n",
            nrow(z_ret), paste(idx_names, collapse = ", "), ncol(z_ret)))

st_eu <- pseudo_normal_block_state(
  z_ret, fix_mean = TRUE, positive = TRUE,
  nested = FALSE, xdelta = 1e-4, tdelta = 1e-4
)
rec_eu <- reconstruct_dynamics(
  st_eu, niter = 25L, sbm = TRUE,
  xdist = FALSE, tdist = FALSE, xl1 = 1.0, tl1 = 0.1
)
print(rec_eu)

prec_eu <- dynamics_precision_matrix(rec_eu)
pcor_eu <- prec_eu$partial_correlation
dimnames(pcor_eu) <- list(idx_names, idx_names)
emp_cor <- stats::cor(log_ret)
cat("\nEmpirical Pearson correlations:\n")
print(round(emp_cor, 3))
cat("\nInferred sparse partial correlations (PseudoNormalBlockState):\n")
print(round(pcor_eu, 3))
cat(sprintf("Minimum eigenvalue of reconstructed precision matrix: %.4f (SPD)\n\n",
            min(eigen(prec_eu$precision, symmetric = TRUE, only.values = TRUE)$values)))

cat("=== 2. Multi-sector Gaussian network (15 assets, 3 sectors, L1 vs SBM) ===\n")
# 3 sectors of 5 assets with sparse intra-sector and bridge precision links
N <- 15L
g_sec <- graph_create(N, directed = FALSE)
e_sec <- rbind(
  cbind(1:4, 2:5), c(5L, 1L), c(1L, 3L),
  cbind(6:9, 7:10), c(10L, 6L), c(6L, 8L),
  cbind(11:14, 12:15), c(15L, 11L), c(11L, 13L),
  c(1L, 6L), c(6L, 11L)
)
graph_add_edges(g_sec, e_sec)
s_sec <- simulate_normal_dynamics(
  g_sec, weights = 0.85, sigma = 0.6, mu = 0,
  n_samples = 800L, burnin = 200L, thin = 3L, mode = "equilibrium"
)
# Standardize as in empirical preprocessing
s_sec_z <- t(scale(t(s_sec)))

st_sec <- pseudo_normal_block_state(
  s_sec_z, fix_mean = TRUE, positive = TRUE,
  nested = TRUE, xdelta = 1e-4, tdelta = 1e-4
)
rec_sec <- reconstruct_dynamics(
  st_sec, niter = 20L, sbm = TRUE,
  xdist = FALSE, tdist = FALSE, xl1 = 1.0, tl1 = 0.1
)
print(rec_sec)

path_sec <- reconstruct_l1_path(
  st_sec, xl1_values = c(1, 5, 15, 30),
  sbm = FALSE, niter = 15L, truth = g_sec
)
print(path_sec)
