test_that("parallel build_candidates is deterministic across thread counts and reconstruct_dynamics returns bs_samples", {
  set.seed(99)
  N <- 50L
  M <- 80L
  g <- graph_create(N, directed = FALSE)
  e1 <- t(combn(1:10, 2))
  e2 <- t(combn(11:20, 2))
  graph_add_edges(g, rbind(e1, e2, c(1L, 11L)))

  S <- simulate_normal_dynamics(g, weights = 0.4, sigma = 0.8, n_samples = M, burnin = 100L, mode = "equilibrium")
  S_z <- t(scale(t(S)))
  rownames(S_z) <- paste0("G", seq_len(N))

  # Verify parallel build_candidates produces identical results across 1 vs 4 threads
  rgraphtool_set_num_threads(1L)
  set.seed(5)
  st1 <- pseudo_normal_block_state(S_z, positive = TRUE, nested = TRUE)
  rec1 <- reconstruct_dynamics(st1, niter = 10L, nsamples = 6L, burnin = 5L, thin = 1L, k_cand = 12L, sbm = TRUE)

  rgraphtool_set_num_threads(4L)
  on.exit(rgraphtool_set_num_threads(1L), add = TRUE)
  set.seed(5)
  st4 <- pseudo_normal_block_state(S_z, positive = TRUE, nested = TRUE)
  rec4 <- reconstruct_dynamics(st4, niter = 10L, nsamples = 6L, burnin = 5L, thin = 1L, k_cand = 12L, sbm = TRUE)

  expect_equal(rec1$entropy, rec4$entropy, tolerance = 1e-10)
  expect_equal(nrow(rec1$edges), nrow(rec4$edges))

  expect_type(rec4$bs_samples, "list")
  expect_length(rec4$bs_samples, 6L)
  mc <- mode_cluster_state(rec4$bs_samples, nested = TRUE)
  expect_gte(mc$K, 1L)
})

test_that("Gram-matrix PseudoNormal matches sample-level PseudoNormal to machine precision", {
  set.seed(2026)
  N <- 20L
  M <- 50L
  g_true <- graph_create(N, directed = FALSE)
  graph_add_edges(g_true, rbind(cbind(1:9, 2:10), c(10L, 1L), cbind(11:19, 12:20), c(20L, 11L)))
  S <- simulate_normal_dynamics(g_true, weights = 0.6, sigma = 0.7, mu = 0.3,
                                n_samples = M, burnin = 100L, mode = "equilibrium")

  # Unshrunk empirical MLE covariance (1/M denominator) + sample means
  mu_hat <- rowMeans(S)
  Sc <- S - mu_hat
  cov_emp <- tcrossprod(Sc) / M

  # Evaluate on an initial random graph with non-zero edges and parameters
  g_init <- random_graph_erdos_renyi(N, p = 0.25, directed = FALSE)
  w_init <- rnorm(graph_num_edges(g_init), 0, 0.25)
  th_init <- cbind(rnorm(N, -0.3, 0.1), rnorm(N, 0.1, 0.1))

  st_samp <- pseudo_normal_block_state(S, graph = g_init, x = w_init, theta = th_init,
                                       positive = TRUE, nested = TRUE,
                                       xdelta = 1e-3, tdelta = c(1e-3, 1e-3))
  st_gram <- pseudo_normal_block_state(cov = cov_emp, n_eff = M, mean = mu_hat,
                                       graph = g_init, x = w_init, theta = th_init,
                                       positive = TRUE, nested = TRUE,
                                       xdelta = 1e-3, tdelta = c(1e-3, 1e-3))

  dl_samp <- dynamics_description_length(st_samp, breakdown = TRUE)
  dl_gram <- dynamics_description_length(st_gram, breakdown = TRUE)
  expect_equal(dl_gram$nll, dl_samp$nll, tolerance = 1e-11)
  expect_equal(dl_gram$total, dl_samp$total, tolerance = 1e-11)

  # Check incremental edge moves match full recomputation
  chk <- rgraphtool:::cpp_dynamics_check_moves(
    st_gram$model_code, st_gram$n_vertices, st_gram$directed, st_gram$self_loops,
    st_gram$s_list, st_gram$f_mat, st_gram$has_zero, st_gram$exposed,
    st_gram$positive, st_gram$pslack, st_gram$x_range[1L], st_gram$x_range[2L],
    st_gram$tmin, st_gram$tmax, st_gram$edges, st_gram$x, st_gram$theta, st_gram$bs,
    st_gram$degree_corrected, st_gram$dist, st_gram$default_prior, 20L, st_gram$q_max
  )
  expect_lt(max(abs(chk)), 1e-11)

  # Reconstruction parity
  set.seed(11)
  rec_s <- reconstruct_dynamics(st_samp, niter = 15L, sbm = TRUE, xdist = FALSE, tdist = FALSE, xl1 = 0.5, tl1 = 0.1)
  set.seed(11)
  rec_g <- reconstruct_dynamics(st_gram, niter = 15L, sbm = TRUE, xdist = FALSE, tdist = FALSE, xl1 = 0.5, tl1 = 0.1)
  expect_equal(rec_g$nll, rec_s$nll, tolerance = 1e-10)
  expect_equal(rec_g$entropy, rec_s$entropy, tolerance = 1e-10)
  expect_equal(nrow(rec_g$edges), nrow(rec_s$edges))
})

test_that("shrinkage_covariance + Gram PseudoNormal work in N >> M regime", {
  set.seed(77)
  N <- 60L
  M <- 18L
  g_true <- graph_create(N, directed = FALSE)
  graph_add_edges(g_true, rbind(cbind(1:19, 2:20), cbind(21:39, 22:40)))
  S <- simulate_normal_dynamics(g_true, weights = 0.5, sigma = 0.7,
                                n_samples = M, burnin = 80L, mode = "equilibrium")
  rownames(S) <- paste0("Gene", seq_len(N))

  sh <- shrinkage_covariance(S, target = "diagonal")
  expect_equal(dim(sh$cov), c(N, N))
  expect_true(isSymmetric(sh$cov, tol = 1e-12))
  expect_true(sh$lambda >= 0 && sh$lambda <= 1)
  eigs <- eigen(sh$cov, symmetric = TRUE, only.values = TRUE)$values
  expect_gt(min(eigs), 0)

  st <- pseudo_normal_block_state(cov = sh$cov, n_eff = sh$n_obs,
                                  positive = TRUE, nested = TRUE,
                                  xdelta = 1e-3, tdelta = c(1e-3, 1e-3))
  rec <- reconstruct_dynamics(st, niter = 15L, k_cand = 20L, sbm = TRUE,
                              xdist = FALSE, tdist = FALSE, xl1 = 0.25, tl1 = 0.1)
  prec <- dynamics_precision_matrix(rec)
  expect_gt(min(eigen(prec$precision, symmetric = TRUE, only.values = TRUE)$values), 0)
})

test_that("M3 active-set refresh recovers suppressor edges missed by marginal screening (F1 > 0.95)", {
  set.seed(314)
  n_tri <- 25L
  N <- 3L * n_tri  # 75 vertices (> 45 -> initial build_candidates screening active)
  Omega <- diag(1.0, N)
  edges_true <- matrix(integer(0), nrow = 0L, ncol = 2L)
  supp_edges <- matrix(integer(0), nrow = n_tri, ncol = 2L)

  for (t in seq_len(n_tri)) {
    a <- 3L * t - 2L
    b <- 3L * t - 1L
    z <- 3L * t
    # Exact suppressor motif: Omega_az * Omega_bz / Omega_zz == Omega_ab => Sigma_ab == 0!
    Omega[a, z] <- Omega[z, a] <- 0.45
    Omega[b, z] <- Omega[z, b] <- -0.45
    Omega[a, b] <- Omega[b, a] <- -0.2025
    edges_true <- rbind(edges_true, c(a, b), c(a, z), c(b, z))
    supp_edges[t, ] <- c(a, b)
  }

  Sigma <- solve(Omega)
  M <- 1500L
  R_chol <- chol(Sigma)
  Z_mat <- matrix(rnorm(M * N), nrow = M, ncol = N)
  S <- t(Z_mat %*% R_chol)

  # Confirm marginal correlation on suppressor pairs is near zero
  cors <- stats::cor(t(S))
  supp_cors <- abs(cors[cbind(supp_edges[, 1L], supp_edges[, 2L])])
  expect_lt(max(supp_cors), 0.12)

  # 1. Fixed initial marginal candidates (refresh_every = 0L) miss suppressor pairs
  st0 <- pseudo_normal_block_state(S, positive = TRUE, nested = TRUE,
                                   xdelta = 1e-3, tdelta = c(1e-3, 1e-3))
  rec0 <- reconstruct_dynamics(st0, niter = 10L, k_cand = 3L, sbm = TRUE,
                               xdist = TRUE, tdist = TRUE, refresh_every = 0L)
  ek <- function(u, v) paste0(pmin(u, v), ":", pmax(u, v))
  supp_keys <- ek(supp_edges[, 1L], supp_edges[, 2L])
  rec0_keys <- ek(rec0$edges$from, rec0$edges$to)
  expect_lt(sum(supp_keys %in% rec0_keys), n_tri / 2)

  # 2. M3 active-set refresh + joint SBM recovers both suppressor pairs and overall graph (F1 > 0.95)
  st_m3 <- pseudo_normal_block_state(S, positive = TRUE, nested = TRUE,
                                     xdelta = 1e-3, tdelta = c(1e-3, 1e-3))
  rec_m3 <- reconstruct_dynamics(st_m3, niter = 10L, k_cand = 3L, sbm = TRUE,
                                 xdist = TRUE, tdist = TRUE, refresh_every = 2L)
  true_keys <- ek(edges_true[, 1L], edges_true[, 2L])
  rec_keys <- ek(rec_m3$edges$from, rec_m3$edges$to)
  tp <- sum(true_keys %in% rec_keys)
  prec <- tp / length(rec_keys)
  recl <- tp / length(true_keys)
  f1 <- 2 * prec * recl / (prec + recl)
  expect_gt(f1, 0.95)
  expect_gte(sum(supp_keys %in% rec_keys), n_tri - 1L)

  # 3. Verify explicit sample-level (non-Gram, M = 60 < 64) active-set refresh codepath
  S_small <- t(matrix(rnorm(60L * N), 60L, N) %*% R_chol)
  st_s <- pseudo_normal_block_state(S_small, positive = TRUE, nested = FALSE)
  rec_s <- reconstruct_dynamics(st_s, niter = 6L, k_cand = 4L, sbm = FALSE,
                                xdist = FALSE, tdist = FALSE, xl1 = 0.5, tl1 = 0.1,
                                refresh_every = 2L)
  expect_true(is.finite(rec_s$entropy))
  expect_equal(rec_s$entropy, dynamics_description_length(rec_s$state), tolerance = 1e-8)
})

test_that("M3 DynBlockPrior merge_split escapes single-block trap and preserves exact DL", {
  set.seed(4242)
  N <- 80L
  Omega <- matrix(0, N, N)
  for (b in 0:3) {
    idx <- (b * 20L + 1L):((b + 1L) * 20L)
    # 5-regular modular blocks with strong weights (-0.25) so MDL recovers >200 intra edges
    for (i in seq_along(idx)) {
      u <- idx[i]
      for (step in c(1L, 2L, 10L)) {
        v <- idx[((i + step - 1L) %% 20L) + 1L]
        Omega[u, v] <- Omega[v, u] <- -0.25
      }
    }
  }
  # Weak inter-block bridges
  for (br in list(c(1L, 21L), c(21L, 41L), c(41L, 61L))) {
    Omega[br[1L], br[2L]] <- Omega[br[2L], br[1L]] <- -0.05
  }
  diag(Omega) <- rowSums(abs(Omega)) + 0.12
  S <- t(matrix(rnorm(600L * N), nrow = 600L, ncol = N) %*% chol(solve(Omega)))
  S_z <- t(scale(t(S)))

  st_no_ms <- pseudo_normal_block_state(S_z, positive = TRUE, nested = TRUE,
                                        blocks = rep(1L, N),
                                        xdelta = 1e-3, tdelta = c(1e-3, 1e-3))
  rec_no_ms <- reconstruct_dynamics(st_no_ms, niter = 10L, k_cand = 20L, sbm = TRUE,
                                    xdist = TRUE, tdist = TRUE, xl1 = 1.0, tl1 = c(0.1, 0.1),
                                    refresh_every = 3L, merge_split = FALSE)

  st_ms <- pseudo_normal_block_state(S_z, positive = TRUE, nested = TRUE,
                                     blocks = rep(1L, N),
                                     xdelta = 1e-3, tdelta = c(1e-3, 1e-3))
  rec_ms <- reconstruct_dynamics(st_ms, niter = 10L, k_cand = 20L, sbm = TRUE,
                                 xdist = TRUE, tdist = TRUE, xl1 = 1.0, tl1 = c(0.1, 0.1),
                                 refresh_every = 3L, merge_split = TRUE)

  B_ms <- length(unique(rec_ms$state$bs[[1L]]))
  expect_gte(B_ms, 3L)
  expect_lt(rec_ms$entropy, rec_no_ms$entropy)
  expect_equal(rec_ms$entropy, dynamics_description_length(rec_ms$state), tolerance = 1e-9)
})

test_that("M3 active-set refresh + merge-split is deterministic across thread counts", {
  set.seed(555)
  N <- 60L
  M <- 100L
  Omega <- diag(1.0, N)
  for (b in 0:1) {
    idx <- (b * 15L + 1L):((b + 1L) * 15L)
    pairs <- t(combn(idx, 2L))
    for (r in seq_len(20L)) {
      u <- pairs[r, 1L]; v <- pairs[r, 2L]
      Omega[u, v] <- Omega[v, u] <- -0.15
    }
  }
  Omega[1L, 16L] <- Omega[16L, 1L] <- -0.15
  S <- t(matrix(rnorm(M * N), nrow = M, ncol = N) %*% chol(solve(Omega)))

  rgraphtool_set_num_threads(1L)
  set.seed(123)
  st1 <- pseudo_normal_block_state(S, positive = TRUE, nested = TRUE, blocks = rep(1L, N))
  rec1 <- reconstruct_dynamics(st1, niter = 8L, k_cand = 10L, sbm = TRUE,
                               xdist = FALSE, tdist = FALSE, xl1 = 0.5, tl1 = 0.1,
                               refresh_every = 2L, merge_split = TRUE)

  rgraphtool_set_num_threads(4L)
  on.exit(rgraphtool_set_num_threads(1L), add = TRUE)
  set.seed(123)
  st4 <- pseudo_normal_block_state(S, positive = TRUE, nested = TRUE, blocks = rep(1L, N))
  rec4 <- reconstruct_dynamics(st4, niter = 8L, k_cand = 10L, sbm = TRUE,
                               xdist = FALSE, tdist = FALSE, xl1 = 0.5, tl1 = 0.1,
                               refresh_every = 2L, merge_split = TRUE)

  expect_equal(rec1$entropy, rec4$entropy, tolerance = 1e-10)
  expect_equal(rec1$nll, rec4$nll, tolerance = 1e-10)
  expect_equal(nrow(rec1$edges), nrow(rec4$edges))
})

test_that("QuantHist incremental MDL matches full recomputation across xdist + tdist", {
  set.seed(808)
  N <- 20L
  M <- 64L
  Omega <- diag(1.0, N)
  for (i in seq(1L, 19L, by = 2L)) {
    Omega[i, i + 1L] <- Omega[i + 1L, i] <- -0.35
  }
  S <- t(matrix(rnorm(M * N), nrow = M, ncol = N) %*% chol(solve(Omega)))
  cov_hat <- stats::cov(t(S)) * (M - 1) / M
  mean_hat <- rowMeans(S)

  st_samp <- pseudo_normal_block_state(
    S, positive = TRUE, nested = TRUE,
    xdelta = 1e-3, tdelta = c(1e-3, 1e-3)
  )
  # Force sample-level evaluation on M=64 by comparing sample vs explicit Gram cov input
  st_gram <- pseudo_normal_block_state(
    cov = cov_hat, n_eff = M, mean = mean_hat,
    positive = TRUE, nested = TRUE,
    xdelta = 1e-3, tdelta = c(1e-3, 1e-3)
  )

  set.seed(42)
  rec_gram <- reconstruct_dynamics(
    st_gram, niter = 10L, k_cand = 8L, sbm = TRUE,
    xdist = TRUE, tdist = TRUE, xl1 = 1.0, tl1 = c(0.1, 0.1),
    refresh_every = 2L, merge_split = TRUE
  )
  expect_equal(rec_gram$entropy, dynamics_description_length(rec_gram$state), tolerance = 1e-9)

  # Exact entropy consistency check on sample state too
  set.seed(42)
  rec_samp <- reconstruct_dynamics(
    st_samp, niter = 10L, k_cand = 8L, sbm = TRUE,
    xdist = TRUE, tdist = TRUE, xl1 = 1.0, tl1 = c(0.1, 0.1),
    refresh_every = 2L, merge_split = TRUE
  )
  expect_equal(rec_samp$entropy, dynamics_description_length(rec_samp$state), tolerance = 1e-9)
  expect_equal(rec_samp$entropy, rec_gram$entropy, tolerance = 1e-9)
  expect_equal(rec_samp$nll, rec_gram$nll, tolerance = 1e-9)
})

test_that("shrinkage_covariance(center = FALSE) shrinks the zero-mean second moment", {
  set.seed(7)
  x <- matrix(rnorm(8 * 30), 8, 30)
  sh <- shrinkage_covariance(x, target = "diagonal", center = FALSE)
  S <- tcrossprod(x) / 30
  expect_equal(diag(sh$cov), diag(S), ignore_attr = TRUE)
  off <- row(S) != col(S)
  expect_equal(sh$cov[off], (1 - sh$lambda) * S[off], ignore_attr = TRUE)
  # invariant to an orthogonal rotation of the observations? second moment yes
  R <- qr.Q(qr(matrix(rnorm(30 * 30), 30)))
  sh_r <- shrinkage_covariance(x %*% R, target = "diagonal", center = FALSE)
  expect_equal(diag(sh_r$cov), diag(sh$cov), ignore_attr = TRUE)
  # the default centred estimator is unchanged
  expect_equal(diag(shrinkage_covariance(x)$cov), apply(x, 1, var), ignore_attr = TRUE)
})
