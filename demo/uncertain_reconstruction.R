# Cookbook demo: Uncertain network reconstruction, latent Poisson multigraphs,
# and Bayesian missing-edge prediction (Peixoto 2018 PRX, Peixoto 2020 PRE).

library(rgraphtool)
set.seed(42)

cat("=== 1. Measured network reconstruction (unknown homogeneous errors) ===\n")
N <- 45L
b_true <- rep(1:3, each = 15L)
all_pairs <- t(utils::combn(N, 2L))
same_blk <- b_true[all_pairs[, 1L]] == b_true[all_pairs[, 2L]]
true_A <- as.integer(stats::runif(nrow(all_pairs)) < ifelse(same_blk, 0.65, 0.02))

# Repeated noisy observations (n = 3 trials, true FN p = 0.20, FP q = 0.05)
n_trials <- 3L
p_fn <- 0.20
q_fp <- 0.05
x_obs <- stats::rbinom(nrow(all_pairs), n_trials, ifelse(true_A == 1L, 1 - p_fn, q_fp))
keep <- x_obs > 0L

g_obs <- graph_create(N, directed = FALSE)
graph_add_edges(
  g_obs,
  all_pairs[keep, , drop = FALSE],
  properties = list(n = rep(n_trials, sum(keep)), x = x_obs[keep])
)

st_meas <- measured_block_state(
  g_obs, n = "n", x = "x",
  n_default = n_trials, x_default = 0L,
  nested = TRUE
)
rec_meas <- reconstruct_network(st_meas, nsamples = 200L, burnin = 150L, thin = 2L)
print(rec_meas)

ek <- function(u, v) paste0(pmin(u, v), ":", pmax(u, v))
idx <- match(ek(all_pairs[, 1L], all_pairs[, 2L]), ek(rec_meas$edges$from, rec_meas$edges$to))
p_hat <- ifelse(is.na(idx), 0, rec_meas$edges$prob[idx])
ev_rec <- evaluate_edge_prediction(p_hat, true_A)
ev_raw <- evaluate_edge_prediction(x_obs / n_trials, true_A)
cat(sprintf(
  "Repeated-trial reconstruction ROC-AUC = %.4f (raw x/n = %.4f) | PR-AUC = %.4f\n",
  ev_rec$roc_auc, ev_raw$roc_auc, ev_rec$pr_auc
))
cat(sprintf(
  "Posterior mean error rates: FN p = %.3f (true %.2f), FP q = %.3f (true %.2f)\n\n",
  mean(rec_meas$p_samples), p_fn, mean(rec_meas$q_samples), q_fp
))

cat("=== 1b. Single-snapshot reconstruction (n = 1, known error rates) ===\n")
x_single <- stats::rbinom(nrow(all_pairs), 1L, ifelse(true_A == 1L, 0.75, 0.05))
keep1 <- x_single > 0L
g_single <- graph_create(N, directed = FALSE)
graph_add_edges(
  g_single,
  all_pairs[keep1, , drop = FALSE],
  properties = list(n = rep(1L, sum(keep1)), x = x_single[keep1])
)
st_fix <- measured_block_state(
  g_single, n = "n", x = "x", n_default = 1L, x_default = 0L,
  lp = log(0.25), lq = log(0.05), nested = TRUE
)
rec_fix <- reconstruct_network(st_fix, nsamples = 200L, burnin = 150L, thin = 2L)
idx1 <- match(ek(all_pairs[, 1L], all_pairs[, 2L]), ek(rec_fix$edges$from, rec_fix$edges$to))
p_hat1 <- ifelse(is.na(idx1), 0, rec_fix$edges$prob[idx1])
ev_fix <- evaluate_edge_prediction(p_hat1, true_A)
ev_raw1 <- evaluate_edge_prediction(x_single, true_A)
cat(sprintf(
  "Single-snapshot ROC-AUC = %.4f vs raw = %.4f | PR-AUC = %.4f vs raw = %.4f\n\n",
  ev_fix$roc_auc, ev_raw1$roc_auc, ev_fix$pr_auc, ev_raw1$pr_auc
))

cat("=== 2. Heterogeneous error reconstruction (MixedMeasuredBlockState) ===\n")
st_mix <- mixed_measured_block_state(
  g_obs, n = "n", x = "x",
  n_default = n_trials, x_default = 0L,
  fn_params = list(alpha = 1, beta = 4),
  fp_params = list(mu = 1, nu = 10),
  nested = TRUE
)
rec_mix <- reconstruct_network(st_mix, nsamples = 100L, burnin = 80L, thin = 2L)
print(rec_mix)
cat(sprintf(
  "Final Beta hyperparameters (alpha, beta, mu, nu): %s\n\n",
  paste(sprintf("%.2f", rec_mix$final_hparams), collapse = ", ")
))

cat("=== 3. Extraneous edge probability estimates (UncertainBlockState) ===\n")
q_est <- pmin(0.98, pmax(0.02, (x_obs[keep] + 0.5) / (n_trials + 1)))
g_unc <- graph_create(N, directed = FALSE)
graph_add_edges(g_unc, all_pairs[keep, , drop = FALSE], properties = list(q = q_est))
st_unc <- uncertain_block_state(g_unc, q = "q", q_default = 0.02, nested = TRUE)
rec_unc <- reconstruct_network(st_unc, nsamples = 100L, burnin = 80L, thin = 2L)
print(rec_unc)
cat("\n")

cat("=== 4. Fast EM latent Poisson multigraph (Peixoto 2020) ===\n")
g_true <- graph_create(N, directed = FALSE)
graph_add_edges(g_true, all_pairs[true_A == 1L, , drop = FALSE])
em_res <- latent_multigraph(g_true, epsilon = 1e-10)
print(em_res)
cat("Top expected edge multiplicities:\n")
print(utils::head(em_res$edges[order(-em_res$edges$weight), ], 5L), row.names = FALSE)
cat("\n")

cat("=== 5. Missing-edge prediction and ROC/PR evaluation ===\n")
true_edges <- all_pairs[true_A == 1L, , drop = FALSE]
non_edges <- all_pairs[true_A == 0L, , drop = FALSE]
hold_n <- min(25L, floor(nrow(true_edges) * 0.25))
hold_idx <- sample(nrow(true_edges), hold_n)
g_train <- graph_create(N, directed = FALSE)
graph_add_edges(g_train, true_edges[-hold_idx, , drop = FALSE])

cands <- rbind(
  true_edges[hold_idx, , drop = FALSE],
  non_edges[sample(nrow(non_edges), hold_n), , drop = FALSE]
)
labels <- c(rep(1L, hold_n), rep(0L, hold_n))
pred <- predict_edges(g_train, cands, nsamples = 60L, burnin = 40L, nested = TRUE)
ev_pred <- evaluate_edge_prediction(pred, labels)
print(pred)
print(ev_pred)

if (interactive() || !identical(Sys.getenv("RGRAPHTOOL_DEMO_PLOT"), "")) {
  op <- graphics::par(mfrow = c(1L, 2L))
  on.exit(graphics::par(op), add = TRUE)
  graphics::plot(
    ev_rec$roc$fpr, ev_rec$roc$tpr, type = "l", lwd = 2, col = "#2b5c8f",
    xlab = "False positive rate", ylab = "True positive rate",
    main = sprintf("Reconstruction ROC (AUC = %.3f)", ev_rec$roc_auc)
  )
  graphics::lines(ev_raw$roc$fpr, ev_raw$roc$tpr, lty = 2, col = "#888888")
  graphics::abline(0, 1, lty = 3, col = "gray60")
  graphics::legend(
    "bottomright", legend = c("SBM reconstruction", "Raw x/n"),
    col = c("#2b5c8f", "#888888"), lty = c(1, 2), lwd = c(2, 1), bty = "n"
  )

  graphics::plot(
    ev_pred$roc$fpr, ev_pred$roc$tpr, type = "l", lwd = 2, col = "#c0392b",
    xlab = "False positive rate", ylab = "True positive rate",
    main = sprintf("Missing-Edge ROC (AUC = %.3f)", ev_pred$roc_auc)
  )
  graphics::abline(0, 1, lty = 3, col = "gray60")
}
