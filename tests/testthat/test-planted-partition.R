test_that("PPBlockState matches pinned graph-tool 3.9 values across all 16 configurations", {
  edges <- rbind(
    c(1L, 2L), c(1L, 2L), c(1L, 3L), c(1L, 1L), c(2L, 3L),
    c(3L, 4L), c(4L, 5L), c(4L, 5L), c(5L, 6L), c(6L, 4L), c(5L, 5L)
  )
  b <- c(1L, 1L, 1L, 2L, 2L, 2L, 3L)

  ref <- list(
    list(dir = FALSE, unif = FALSE, dc = TRUE,  dk = "distributed", S = 30.045371343521),
    list(dir = FALSE, unif = FALSE, dc = TRUE,  dk = "uniform",     S = 30.323239794524),
    list(dir = FALSE, unif = FALSE, dc = FALSE, dk = "distributed", S = 27.804380235024),
    list(dir = FALSE, unif = FALSE, dc = FALSE, dk = "uniform",     S = 27.804380235024),
    list(dir = FALSE, unif = TRUE,  dc = TRUE,  dk = "distributed", S = 31.312410400664),
    list(dir = FALSE, unif = TRUE,  dc = TRUE,  dk = "uniform",     S = 31.590278851668),
    list(dir = FALSE, unif = TRUE,  dc = FALSE, dk = "distributed", S = 29.071419292167),
    list(dir = FALSE, unif = TRUE,  dc = FALSE, dk = "uniform",     S = 29.071419292167),
    list(dir = TRUE,  unif = FALSE, dc = TRUE,  dk = "distributed", S = 36.369463674785),
    list(dir = TRUE,  unif = FALSE, dc = TRUE,  dk = "uniform",     S = 36.924625112371),
    list(dir = TRUE,  unif = FALSE, dc = FALSE, dk = "distributed", S = 34.042704860064),
    list(dir = TRUE,  unif = FALSE, dc = FALSE, dk = "uniform",     S = 34.042704860064),
    list(dir = TRUE,  unif = TRUE,  dc = TRUE,  dk = "distributed", S = 37.636502731928),
    list(dir = TRUE,  unif = TRUE,  dc = TRUE,  dk = "uniform",     S = 38.191664169514),
    list(dir = TRUE,  unif = TRUE,  dc = FALSE, dk = "distributed", S = 35.309743917207),
    list(dir = TRUE,  unif = TRUE,  dc = FALSE, dk = "uniform",     S = 35.309743917207)
  )

  g_u <- graph_add_edges(graph_create(7, directed = FALSE), edges)
  g_d <- graph_add_edges(graph_create(7, directed = TRUE), edges)

  for (r in ref) {
    g <- if (r$dir) g_d else g_u
    st <- pp_block_state(g, blocks = b, uniform = r$unif,
                         degree_corrected = r$dc, degree_prior = r$dk)
    expect_equal(pp_description_length(st), r$S, tolerance = 1e-10)
    set.seed(11L)
    errs <- rgraphtool:::cpp_pp_check_moves(
      rgraphtool:::graph_pointer(g), st$blocks, st$uniform, st$degree_corrected,
      r$dk == "distributed", 30L, 1000
    )
    expect_lt(max(errs), 1e-9)
  }
})

test_that("planted-partition recovers assortative communities and rejects disassortative ones", {
  set.seed(42L)
  # 1. Strongly assortative 3-community graph
  n <- 36L
  truth <- rep(1:3, each = 12L)
  ed_a <- list()
  for (i in 1:(n - 1L)) {
    for (j in (i + 1L):n) {
      p <- if (truth[i] == truth[j]) 0.70 else 0.03
      if (stats::runif(1) < p) ed_a[[length(ed_a) + 1L]] <- c(i, j)
    }
  }
  g_assort <- graph_add_edges(graph_create(n, directed = FALSE), do.call(rbind, ed_a))

  st_pp <- minimize_pp_blockmodel_dl(g_assort)
  expect_equal(st_pp$B, 3L)
  expect_gt(partition_overlap(block_membership(st_pp), truth, norm = TRUE), 0.95)
  expect_equal(sum(block_sizes(st_pp)), n)
  asrt <- pp_assortativity(st_pp)
  expect_gt(asrt$fraction_in, 0.85)
  expect_equal(asrt$E_in + asrt$E_out, asrt$E_total)
  expect_output(print(st_pp), "planted-partition block state")

  # Assortative model comparison favors PP on a pure planted-partition graph
  cmp_a <- compare_assortative_models(g_assort, method = "map")
  expect_s3_class(cmp_a, "rgraphtool_assortative_model_comparison")
  expect_true(cmp_a$best_model %in% c("pp", "pp_uniform"))
  expect_equal(sum(cmp_a$summary$posterior_weight), 1, tolerance = 1e-9)
  expect_output(print(cmp_a), "Best model:")

  # 2. Mixed assortative + bipartite 4-block network:
  #    blocks 1 & 2 are dense assortative cliques, while blocks 3 & 4 form a
  #    dense bipartite pair. No single (E_in, E_out) rate can explain both, so
  #    the general SBM decisively beats PP.
  ed_b <- list()
  truth_b <- rep(1:4, each = 10L)
  P_mat <- rbind(
    c(0.75, 0.02, 0.02, 0.02),
    c(0.02, 0.75, 0.02, 0.02),
    c(0.02, 0.02, 0.02, 0.75),
    c(0.02, 0.02, 0.75, 0.02)
  )
  for (i in 1:39) {
    for (j in (i + 1L):40) {
      if (stats::runif(1) < P_mat[truth_b[i], truth_b[j]]) {
        ed_b[[length(ed_b) + 1L]] <- c(i, j)
      }
    }
  }
  g_nonassort <- graph_add_edges(graph_create(40L, directed = FALSE), do.call(rbind, ed_b))
  cmp_b <- compare_assortative_models(g_nonassort, method = "map")
  expect_true(cmp_b$best_model %in% c("sbm_flat", "sbm_nested"))
})

test_that("PPBlockState supports MCMC, posterior sampling, model_evidence, and validation", {
  set.seed(9L)
  edges <- rbind(t(utils::combn(1:8, 2)), t(utils::combn(9:16, 2)), c(8, 9), c(1, 10))
  g <- graph_add_edges(graph_create(16L, directed = FALSE), edges)
  st <- minimize_pp_blockmodel_dl(g)

  sw <- mcmc_sweep(st, beta = 1.0, niter = 10L)
  expect_s3_class(sw$state, "rgraphtool_pp_block_state")
  expect_equal(sw$entropy, pp_description_length(sw$state), tolerance = 1e-9)

  eq <- mcmc_equilibrate(st, wait = 5L, nbreaks = 1L, max_niter = 15L)
  expect_s3_class(eq$state, "rgraphtool_pp_block_state")

  an <- mcmc_anneal(st, beta_range = c(1, 5), niter = 3L,
                    mcmc_equilibrate_args = list(wait = 3L, max_niter = 6L))
  expect_s3_class(an$state, "rgraphtool_pp_block_state")

  post <- sample_posterior(st, nsamples = 20L, burnin = 15L, thin = 2L)
  expect_equal(length(post$samples), 20L)
  expect_false(post$nested)

  ev <- model_evidence(st, nsamples = 20L, burnin = 15L, thin = 2L)
  expect_s3_class(ev, "rgraphtool_model_evidence")
  expect_true(is.finite(ev$log_evidence))

  cmp_p <- compare_assortative_models(g, models = c("pp", "sbm_flat"),
                                      method = "posterior", nsamples = 20L,
                                      burnin = 15L, thin = 2L)
  expect_true(all(is.finite(cmp_p$summary$neg_log_evidence)))

  # Validation & graph mutation
  g0 <- graph_create(4, directed = FALSE)
  expect_true(is.finite(pp_description_length(pp_block_state(g0))))
  expect_error(pp_block_state(g, uniform = NA), "TRUE or FALSE")
  expect_error(compare_assortative_models(g, models = "invalid"), "models must be")
  graph_add_edges(g, rbind(c(2L, 12L)))
  expect_error(pp_description_length(st), "changed")
})

test_that("weighted PPBlockState matches B=1 WeightedBlockState and exact incremental moves", {
  set.seed(501L)
  edges <- rbind(
    c(1L, 2L), c(1L, 2L), c(1L, 3L), c(1L, 1L), c(2L, 3L),
    c(3L, 4L), c(4L, 5L), c(4L, 5L), c(5L, 6L), c(6L, 4L), c(5L, 5L)
  )
  m <- nrow(edges)
  b1 <- rep(1L, 7L)
  b3 <- c(1L, 1L, 1L, 2L, 2L, 3L, 3L)
  w_norm <- stats::rnorm(m, 0.4, 1.1)
  w_exp <- exp(0.3 * w_norm)
  w_pois <- stats::rpois(m, lambda = 4.0)

  specs <- list(
    list(w = w_norm, tp = "real-normal", pr = "microcanonical"),
    list(w = w_norm, tp = "real-normal", pr = "bayesian"),
    list(w = w_exp,  tp = "real-exponential", pr = "microcanonical"),
    list(w = w_pois, tp = "discrete-poisson", pr = "microcanonical")
  )

  for (dir in c(FALSE, TRUE)) {
    g <- graph_add_edges(graph_create(7L, directed = dir), edges)
    for (sp in specs) {
      for (dc in c(TRUE, FALSE)) {
        ws <- weighted_block_state(g, weights = sp$w, types = sp$tp, params = sp$pr,
                                   blocks = b1, degree_corrected = dc)
        dl_w <- weighted_description_length(ws)
        for (unif in c(FALSE, TRUE)) {
          ps <- pp_block_state(g, blocks = b1, uniform = unif,
                               degree_corrected = dc,
                               weights = sp$w, types = sp$tp, params = sp$pr)
          expect_equal(pp_description_length(ps), dl_w, tolerance = 1e-10)
        }
      }
    }

    w_multi <- cbind(w_norm, w_pois)
    tp_multi <- c("real-normal", "discrete-poisson")
    pr_multi <- list(c(NaN, NaN, NaN, NaN), c(NaN, NaN))
    for (unif in c(FALSE, TRUE)) {
      for (dc in c(TRUE, FALSE)) {
        err <- rgraphtool:::cpp_pp_check_moves(
          rgraphtool:::graph_pointer(g), b3, unif, dc, TRUE, 40L, 1000,
          w_multi, tp_multi, pr_multi
        )
        expect_lt(max(err), 1e-8)
      }
    }
  }
})

test_that("weighted PP recovers dense weight-only communities and supports compare_assortative_models", {
  set.seed(502L)
  # 1. Uniform dense topology (p = 0.6 everywhere): signal is purely in edge weights
  n <- 30L
  truth <- rep(1:3, each = 10L)
  ed <- list()
  w <- numeric(0)
  for (i in 1:(n - 1L)) {
    for (j in (i + 1L):n) {
      if (stats::runif(1) < 0.6) {
        ed[[length(ed) + 1L]] <- c(i, j)
        same <- (truth[i] == truth[j])
        w <- c(w, stats::rnorm(1, mean = if (same) 2.2 else -0.6, sd = 0.45))
      }
    }
  }
  g_dense <- graph_add_edges(
    graph_create(n, directed = FALSE),
    do.call(rbind, ed),
    properties = list(w = w)
  )

  st_w <- minimize_pp_blockmodel_dl(g_dense, weights = "w", types = "real-normal")
  expect_true(isTRUE(st_w$weighted))
  expect_equal(st_w$B, 3L)
  expect_gt(partition_overlap(block_membership(st_w), truth, norm = TRUE), 0.90)
  asrt <- pp_assortativity(st_w)
  expect_gt(asrt$weight_in_mean, 1.5)
  expect_lt(asrt$weight_out_mean, 0.0)
  expect_output(print(st_w), "weighted \\(real-normal\\)")

  cmp_dense <- compare_assortative_models(
    g_dense, weights = "w", types = "real-normal",
    method = "posterior", nsamples = 25L, burnin = 15L, thin = 2L
  )
  expect_true(cmp_dense$best_model %in% c("pp", "pp_uniform"))
  expect_output(print(cmp_dense), "weighted \\(real-normal\\)")

  # MCMC / posterior / model_evidence preserve weights
  sw <- mcmc_sweep(st_w, beta = 1.0, niter = 8L)
  expect_true(isTRUE(sw$state$weighted))
  expect_equal(sw$entropy, pp_description_length(sw$state), tolerance = 1e-8)
  post <- sample_posterior(st_w, nsamples = 15L, burnin = 10L, thin = 2L)
  expect_equal(length(post$samples), 15L)
  ev <- model_evidence(st_w, nsamples = 15L, burnin = 10L, thin = 2L)
  expect_true(is.finite(ev$log_evidence))

  # 2. Structured inter-module antagonism vs background: general weighted SBM beats PP
  n4 <- 32L
  truth4 <- rep(1:4, each = 8L)
  mu_mat <- rbind(
    c( 2.0, -2.5,  0.0,  0.0),
    c(-2.5,  2.0,  0.0,  0.0),
    c( 0.0,  0.0,  2.0, -2.5),
    c( 0.0,  0.0, -2.5,  2.0)
  )
  ed4 <- list()
  w4 <- numeric(0)
  for (i in 1:(n4 - 1L)) {
    for (j in (i + 1L):n4) {
      if (stats::runif(1) < 0.7) {
        ed4[[length(ed4) + 1L]] <- c(i, j)
        w4 <- c(w4, stats::rnorm(1, mean = mu_mat[truth4[i], truth4[j]], sd = 0.35))
      }
    }
  }
  g_antag <- graph_add_edges(
    graph_create(n4, directed = FALSE),
    do.call(rbind, ed4),
    properties = list(w = w4)
  )
  cmp_antag <- compare_assortative_models(
    g_antag, weights = "w", types = "real-normal", method = "map"
  )
  expect_true(cmp_antag$best_model %in% c("sbm_flat", "sbm_nested"))
})

test_that("layered and layered-weighted PPBlockState match B=1 and C=1 references and incremental moves", {
  set.seed(503L)
  edges <- rbind(
    c(1L, 2L), c(1L, 2L), c(1L, 3L), c(1L, 1L), c(2L, 3L),
    c(3L, 4L), c(4L, 5L), c(4L, 5L), c(5L, 6L), c(6L, 4L), c(5L, 5L)
  )
  m <- nrow(edges)
  # Vertex 7 is isolated; multi-layer parallel + self-loops present
  lay1 <- rep(1L, m)
  lay2 <- rep(1:2, length.out = m)
  w_norm <- stats::rnorm(m, 0.3, 1.2)
  w_exp <- exp(0.25 * w_norm)
  w_pois <- stats::rpois(m, lambda = 4.5)
  b1 <- rep(1L, 7L)
  b3 <- c(1L, 1L, 1L, 2L, 2L, 3L, 3L)

  for (dir in c(FALSE, TRUE)) {
    g <- graph_add_edges(graph_create(7L, directed = dir), edges)

    # C = 1 exact parity with single-layer PPBlockState
    for (unif in c(FALSE, TRUE)) {
      p_single <- pp_block_state(g, blocks = b3, uniform = unif, weights = w_norm, types = "real-normal")
      p_lay1 <- pp_block_state(g, blocks = b3, uniform = unif, weights = w_norm, types = "real-normal",
                               layers = lay1, independent = FALSE)
      expect_equal(pp_description_length(p_lay1), pp_description_length(p_single), tolerance = 1e-10)
    }

    # B = 1 exact parity with LayeredBlockState across unweighted & weighted configs
    specs <- list(
      list(w = NULL, tp = NULL, pr = "microcanonical"),
      list(w = w_norm, tp = "real-normal", pr = "microcanonical"),
      list(w = w_exp,  tp = "real-exponential", pr = "microcanonical"),
      list(w = w_pois, tp = "discrete-poisson", pr = "microcanonical")
    )
    for (sp in specs) {
      for (indep in c(FALSE, TRUE)) {
        for (dc in c(TRUE, FALSE)) {
          ls <- layered_block_state(g, layers = lay2, independent = indep,
                                    blocks = b1, degree_corrected = dc,
                                    weights = sp$w, types = sp$tp, params = sp$pr)
          dl_l <- layered_description_length(ls)
          for (unif in c(FALSE, TRUE)) {
            ps <- pp_block_state(g, blocks = b1, uniform = unif, degree_corrected = dc,
                                 weights = sp$w, types = sp$tp, params = sp$pr,
                                 layers = lay2, independent = indep)
            expect_equal(pp_description_length(ps), dl_l, tolerance = 1e-10)
          }
        }
      }
    }

    # All 8 (unif x indep x dc) incremental move checks with multi-covariate weights
    w_multi <- cbind(w_norm, w_pois)
    tp_multi <- c("real-normal", "discrete-poisson")
    pr_multi <- list(c(NaN, NaN, NaN, NaN), c(NaN, NaN))
    for (unif in c(FALSE, TRUE)) {
      for (indep in c(FALSE, TRUE)) {
        for (dc in c(TRUE, FALSE)) {
          err <- rgraphtool:::cpp_pp_check_moves(
            rgraphtool:::graph_pointer(g), b3, unif, dc, TRUE, 45L, 1000,
            w_multi, tp_multi, pr_multi, lay2, indep
          )
          expect_lt(max(err), 1e-8)
        }
      }
    }
  }
})

test_that("layered PP recovers multi-layer communities and integrates with compare_layered_models", {
  set.seed(504L)
  n <- 24L
  truth <- rep(1:2, each = 12L)
  ed_list <- list()
  lay <- integer(0)
  w_cons <- numeric(0)
  w_rew <- numeric(0)

  for (i in 1:(n - 1L)) {
    for (j in (i + 1L):n) {
      same <- (truth[i] == truth[j])
      if (stats::runif(1) < 0.7) {
        for (l in 1:2) {
          ed_list[[length(ed_list) + 1L]] <- c(i, j)
          lay <- c(lay, l)
          mu_c <- if (same) 2.4 else -1.2
          mu_r <- if (l == 1L) (if (same) 2.4 else -1.2) else (if (same) -2.4 else 1.2)
          w_cons <- c(w_cons, stats::rnorm(1, mu_c, 0.35))
          w_rew <- c(w_rew, stats::rnorm(1, mu_r, 0.35))
        }
      }
    }
  }
  g <- graph_add_edges(
    graph_create(n, directed = FALSE),
    do.call(rbind, ed_list),
    properties = list(layer = lay, w_cons = w_cons, w_rew = w_rew)
  )

  st_lp <- minimize_pp_blockmodel_dl(
    g, weights = "w_cons", types = "real-normal", layers = "layer", independent = FALSE
  )
  expect_true(isTRUE(st_lp$layered))
  expect_true(isTRUE(st_lp$weighted))
  expect_equal(st_lp$C, 2L)
  expect_gt(partition_overlap(block_membership(st_lp), truth, norm = TRUE), 0.90)
  asrt <- pp_assortativity(st_lp)
  expect_equal(dim(asrt$err_layer), c(st_lp$B, 2L))
  expect_equal(dim(asrt$weight_in_mean_layer), c(2L, 1L))
  expect_output(print(st_lp), "2 layers \\(edge-covariates\\)")

  sw <- mcmc_sweep(st_lp, beta = 1.0, niter = 6L)
  expect_true(isTRUE(sw$state$layered))
  expect_true(isTRUE(sw$state$weighted))
  expect_equal(sw$entropy, pp_description_length(sw$state), tolerance = 1e-8)

  cmp_cons <- compare_layered_models(
    g, layers = "layer",
    models = c("edge_covariates", "independent_layers",
               "pp_edge_covariates", "pp_independent_layers"),
    weights = "w_cons", types = "real-normal",
    nested = TRUE, method = "posterior", nsamples = 25L, burnin = 15L, thin = 2L
  )
  expect_equal(cmp_cons$best_model, "pp_edge_covariates")

  cmp_rew <- compare_layered_models(
    g, layers = "layer",
    models = c("pp_edge_covariates", "pp_independent_layers",
               "edge_covariates", "independent_layers"),
    weights = "w_rew", types = "real-normal",
    nested = TRUE, method = "posterior", nsamples = 25L, burnin = 15L, thin = 2L
  )
  expect_true(cmp_rew$best_model %in% c("pp_independent_layers", "independent_layers"))
  expect_true(cmp_rew$summary$independent[1L])
})

test_that("independent-layer PP state on an edgeless graph does not crash", {
  st <- pp_block_state(graph_create(5L), layers = integer(0), independent = TRUE)
  expect_true(is.finite(pp_description_length(st)))
})
