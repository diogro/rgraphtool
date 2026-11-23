# Cookbook demo: Network reconstruction from synthetic dynamics and behavior
# (Ising, epidemics, linear normal SDE, Potts) and L1 regularization comparison
# (Peixoto 2019 PRL, Peixoto 2025 PRX).

library(rgraphtool)
set.seed(2026)

eval_graph <- function(rec, g_true) {
  N <- graph_num_vertices(g_true)
  dir <- isTRUE(g_true$directed)
  ek <- function(u, v) if (!dir) paste0(pmin(u, v), ":", pmax(u, v)) else paste0(u, ":", v)
  te <- graph_edges(g_true)
  tk <- unique(ek(te[, 1L], te[, 2L]))
  if (dir) {
    ap <- expand.grid(u = seq_len(N), v = seq_len(N))
    if (!rec$state$self_loops) ap <- ap[ap$u != ap$v, , drop = FALSE]
    ap <- cbind(ap$u, ap$v)
  } else {
    ap <- t(utils::combn(N, 2L))
  }
  pk <- ek(ap[, 1L], ap[, 2L])
  y <- as.integer(pk %in% tk)
  rk <- ek(rec$edges$from, rec$edges$to)
  idx <- match(pk, rk)
  sc <- ifelse(is.na(idx), 0, abs(rec$edges$weight[idx]) * rec$edges$prob[idx])
  ev <- evaluate_edge_prediction(sc, y, threshold = 1e-6)
  list(roc_auc = ev$roc_auc, pr_auc = ev$pr_auc, f1 = ev$f1,
       n_pred = nrow(rec$edges), n_true = sum(y))
}

# Build a 2-community undirected sparse ground-truth graph (20 vertices, 26 edges)
N <- 20L
g_u <- graph_create(N, directed = FALSE)
e_u <- rbind(
  cbind(1:9, 2:10), c(10L, 1L), c(1L, 5L), c(3L, 8L),
  cbind(11:19, 12:20), c(20L, 11L), c(11L, 15L), c(13L, 18L),
  c(1L, 11L), c(6L, 16L)
)
graph_add_edges(g_u, e_u)

cat("=== 1. Equilibrium Ising reconstruction: SBM+MDL vs. L1 regularization ===\n")
s_pi <- simulate_ising(g_u, weights = 0.55, theta = 0,
                       n_samples = 600L, burnin = 300L, thin = 4L,
                       mode = "equilibrium")
st_pi <- pseudo_ising_block_state(s_pi, nested = TRUE, xdelta = 1e-4, tdelta = 1e-4)
rec_pi <- reconstruct_dynamics(st_pi, niter = 18L, sbm = TRUE,
                               xdist = FALSE, tdist = FALSE, xl1 = 1.0, tl1 = 0.1)
ev_pi <- eval_graph(rec_pi, g_u)
print(rec_pi)
cat(sprintf("SBM+L1 recovery: ROC-AUC=%.4f, PR-AUC=%.4f, F1=%.4f (edges: %d/%d)\n\n",
            ev_pi$roc_auc, ev_pi$pr_auc, ev_pi$f1, ev_pi$n_pred, ev_pi$n_true))

cat("--- Pure L1 regularization path (sbm = FALSE) ---\n")
l1_path <- reconstruct_l1_path(st_pi, xl1_values = c(1, 5, 10, 20),
                               sbm = FALSE, niter = 12L, truth = g_u)
print(l1_path)
cat("\n")

cat("=== 2. Directed epidemic cascade reconstruction (SIS, 20 cascades) ===\n")
g_d <- graph_create(16L, directed = TRUE)
e_d <- rbind(cbind(1:15, 2:16), c(16L, 1L), c(1L, 5L), c(5L, 9L), c(9L, 13L), c(13L, 1L))
graph_add_edges(g_d, e_d)
casc <- simulate_epidemics(g_d, transmission = 0.45, spontaneous = 0.01,
                           recovery = 0.3, model = "SIS",
                           n_steps = 60L, n_cascades = 20L)
st_ep <- epidemics_block_state(casc, directed = TRUE, nested = FALSE,
                               xdelta = 1e-4, tdelta = 1e-4)
rec_ep <- reconstruct_dynamics(st_ep, niter = 15L, sbm = TRUE,
                               xdist = FALSE, tdist = FALSE, xl1 = 1.0, tl1 = 0.1)
ev_ep <- eval_graph(rec_ep, g_d)
tp_ep <- dynamics_transmission_probs(rec_ep)
cat(sprintf(
  "Epidemics recovery: ROC-AUC=%.4f, F1=%.4f (edges: %d/%d) | median t_hat=%.3f (true 0.45)\n\n",
  ev_ep$roc_auc, ev_ep$f1, ev_ep$n_pred, ev_ep$n_true,
  stats::median(tp_ep$edges$transmission_prob)
))

cat("=== 3. Directed linear Gaussian dynamics (LinearNormalBlockState) ===\n")
g_ln <- graph_create(16L, directed = TRUE)
e_ln <- rbind(e_d, cbind(1:16, 1:16))
w_ln <- c(rep(0.35, nrow(e_d)), rep(-0.6, 16L))
graph_add_edges(g_ln, e_ln, properties = list(w = w_ln))
s_ln <- simulate_normal_dynamics(g_ln, weights = "w", sigma = 0.5,
                                 n_samples = 600L, burnin = 100L, thin = 1L,
                                 mode = "linear")
st_ln <- linear_normal_block_state(s_ln, directed = TRUE, self_loops = TRUE,
                                   nested = FALSE, xdelta = 1e-4, tdelta = 1e-4)
rec_ln <- reconstruct_dynamics(st_ln, niter = 15L, sbm = TRUE,
                               xdist = FALSE, tdist = FALSE, xl1 = 1.0, tl1 = 0.1)
ev_ln <- eval_graph(rec_ln, g_ln)
cat(sprintf("LinearNormal recovery: ROC-AUC=%.4f, F1=%.4f (edges: %d/%d)\n\n",
            ev_ln$roc_auc, ev_ln$f1, ev_ln$n_pred, ev_ln$n_true))

cat("=== 4. Categorical Potts / voting behavior (q = 3, PseudoPottsBlockState) ===\n")
s_pt <- simulate_potts(g_u, weights = 1.2, f = diag(3),
                       n_samples = 400L, burnin = 200L, thin = 4L,
                       mode = "equilibrium")
st_pt <- pseudo_potts_block_state(s_pt, f = diag(3), nested = FALSE,
                                  xdelta = 1e-4, tdelta = 1e-4)
rec_pt <- reconstruct_dynamics(st_pt, niter = 12L, sbm = TRUE,
                               xdist = FALSE, tdist = FALSE, xl1 = 1.0, tl1 = 0.1)
ev_pt <- eval_graph(rec_pt, g_u)
cat(sprintf("PseudoPotts recovery: ROC-AUC=%.4f, F1=%.4f (edges: %d/%d)\n",
            ev_pt$roc_auc, ev_pt$f1, ev_pt$n_pred, ev_pt$n_true))
