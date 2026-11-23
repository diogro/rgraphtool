# Tests for Milestone 10: Uncertain network reconstruction, measured networks,
# heterogeneous errors, extraneous probabilities, latent Poisson multigraphs,
# and Bayesian edge prediction.

make_m10_ref_graph <- function(directed = FALSE) {
  edges <- matrix(c(
    1L, 2L,
    2L, 3L,
    3L, 1L,
    4L, 5L,
    5L, 6L,
    6L, 4L,
    3L, 4L
  ), ncol = 2L, byrow = TRUE)
  n_vals <- c(3, 3, 2, 4, 3, 2, 1)
  x_vals <- c(3, 2, 2, 4, 3, 1, 1)
  q_vals <- c(0.95, 0.85, 0.90, 0.99, 0.80, 0.70, 0.40)
  g <- graph_create(7L, directed = directed)
  graph_add_edges(g, edges, properties = list(n = n_vals, x = x_vals, q = q_vals))
  g
}

test_that("measured_block_state matches pinned graph-tool 3.9 to machine precision", {
  b_vals <- c(1L, 1L, 1L, 2L, 2L, 2L, 3L)
  bs_vals <- list(b_vals, c(1L, 1L, 2L), c(1L, 1L))

  # Pinned graph-tool 3.9 reference values
  cases <- list(
    list(dir = FALSE, nest = FALSE, dc = TRUE,  fix = FALSE,
         tot = 34.16353426089133, meas = 8.585449395995033, dens = 6.912318298849235,
         pp = c(3.5, 18.5), qp = c(1.2, 31.0)),
    list(dir = FALSE, nest = FALSE, dc = TRUE,  fix = TRUE,
         tot = 31.513476075032976, meas = 5.935391210136675),
    list(dir = FALSE, nest = FALSE, dc = FALSE, fix = FALSE,
         tot = 32.49378358011404, meas = 8.585449395995033, dens = 6.912318298849235,
         pp = c(3.5, 18.5), qp = c(1.2, 31.0)),
    list(dir = FALSE, nest = TRUE,  dc = TRUE,  fix = FALSE,
         tot = 38.23952968444524, meas = 8.585449395995033, dens = 6.912318298849235,
         pp = c(3.5, 18.5), qp = c(1.2, 31.0),
         ep03 = -5.165322408106305, ep23 = -1.652901361179718),
    list(dir = FALSE, nest = TRUE,  dc = TRUE,  fix = TRUE,
         tot = 35.58947149858688, meas = 5.935391210136675),
    list(dir = TRUE,  nest = FALSE, dc = TRUE,  fix = FALSE,
         tot = 40.11297374961078, meas = 9.611012680325507, dens = 7.937881583179709,
         pp = c(3.5, 18.5), qp = c(1.2, 73.0)),
    list(dir = TRUE,  nest = FALSE, dc = FALSE, fix = TRUE,
         tot = 38.25187257010828, meas = 8.0897095744138),
    list(dir = TRUE,  nest = TRUE,  dc = TRUE,  fix = FALSE,
         tot = 44.50196905360076, meas = 9.611012680325507, dens = 7.937881583179709,
         pp = c(3.5, 18.5), qp = c(1.2, 73.0),
         ep03 = -5.168320151366225, ep23 = -1.1691886326738952)
  )

  q_edges <- matrix(c(1L, 4L, 3L, 4L), ncol = 2L, byrow = TRUE)

  for (cs in cases) {
    g <- make_m10_ref_graph(cs$dir)
    bl <- if (cs$nest) bs_vals else b_vals
    st <- if (!cs$fix) {
      measured_block_state(
        g, n = "n", x = "x", n_default = 2L, x_default = 0L,
        fn_params = list(alpha = 1.5, beta = 2.5),
        fp_params = list(mu = 1.2, nu = 3.0),
        nested = cs$nest, blocks = bl, degree_corrected = cs$dc, self_loops = FALSE
      )
    } else {
      measured_block_state(
        g, n = "n", x = "x", n_default = 2L, x_default = 0L,
        lp = log(0.1), lq = log(0.05),
        nested = cs$nest, blocks = bl, degree_corrected = cs$dc, self_loops = FALSE
      )
    }
    dl <- uncertain_description_length(st, breakdown = TRUE)
    expect_equal(dl$total, cs$tot, tolerance = 1e-11)
    expect_equal(dl$measurement, cs$meas, tolerance = 1e-11)
    expect_equal(uncertain_description_length(st), cs$tot, tolerance = 1e-11)
    if (!is.null(cs$dens)) {
      expect_equal(uncertain_description_length(st, sbm = FALSE, density = TRUE, aE = 2.0),
                   cs$dens, tolerance = 1e-11)
    }
    if (!is.null(cs$pp)) {
      ep <- measured_error_posterior(st)
      expect_equal(c(ep$p_alpha, ep$p_beta), cs$pp, tolerance = 1e-12)
      expect_equal(c(ep$q_alpha, ep$q_beta), cs$qp, tolerance = 1e-12)
    }
    if (!is.null(cs$ep03)) {
      lps <- edge_probabilities(st, q_edges, log = TRUE)
      expect_equal(lps[1L], cs$ep03, tolerance = 1e-11)
      expect_equal(lps[2L], cs$ep23, tolerance = 1e-11)
      expect_equal(edge_probabilities(st, q_edges, log = FALSE), exp(lps), tolerance = 1e-12)
    }
  }
})

test_that("mixed_measured, uncertain, and latent_multigraph_block_state match graph-tool 3.9", {
  b_vals <- c(1L, 1L, 1L, 2L, 2L, 2L, 3L)
  bs_vals <- list(b_vals, c(1L, 1L, 2L), c(1L, 1L))
  q_edges <- matrix(c(1L, 4L, 3L, 4L), ncol = 2L, byrow = TRUE)

  # MixedMeasuredBlockState
  mixed_cases <- list(
    list(dir = FALSE, nest = FALSE, tot = 34.54606084884234, meas = 8.967975983946037),
    list(dir = FALSE, nest = TRUE,  tot = 38.62205627239624, meas = 8.967975983946037,
         ep03 = -5.237708946771244, ep23 = -2.386029349916799),
    list(dir = TRUE,  nest = FALSE, tot = 43.32559557918139, meas = 12.823634509896113),
    list(dir = TRUE,  nest = TRUE,  tot = 47.71459088317137, meas = 12.823634509896113,
         ep03 = -5.19448249990979, ep23 = -2.561542487253084)
  )
  for (cs in mixed_cases) {
    g <- make_m10_ref_graph(cs$dir)
    st <- mixed_measured_block_state(
      g, n = "n", x = "x", n_default = 2L, x_default = 0L,
      fn_params = list(alpha = 1.5, beta = 10.0),
      fp_params = list(mu = 1.2, nu = 12.0),
      nested = cs$nest, blocks = if (cs$nest) bs_vals else b_vals
    )
    dl <- uncertain_description_length(st, breakdown = TRUE)
    expect_equal(dl$total, cs$tot, tolerance = 1e-11)
    expect_equal(dl$measurement, cs$meas, tolerance = 1e-11)
    if (!is.null(cs$ep03)) {
      lps <- edge_probabilities(st, q_edges, log = TRUE)
      expect_equal(lps[1L], cs$ep03, tolerance = 1e-11)
      expect_equal(lps[2L], cs$ep23, tolerance = 1e-11)
    }
  }

  # UncertainBlockState
  unc_cases <- list(
    list(dir = FALSE, nest = FALSE, tot = 14.17377970375718, meas = -11.40430516113912),
    list(dir = FALSE, nest = TRUE,  tot = 18.249775127311086, meas = -11.40430516113912,
         ep23 = -4.0243353357305445),
    list(dir = TRUE,  nest = FALSE, tot = 14.065068311229673, meas = -16.436892758055603),
    list(dir = TRUE,  nest = TRUE,  tot = 18.45406361521965, meas = -16.436892758055603,
         ep23 = -3.4409041758322094)
  )
  for (cs in unc_cases) {
    g <- make_m10_ref_graph(cs$dir)
    st <- uncertain_block_state(
      g, q = "q", q_default = 0.02,
      nested = cs$nest, blocks = if (cs$nest) bs_vals else b_vals
    )
    dl <- uncertain_description_length(st, breakdown = TRUE)
    expect_equal(dl$total, cs$tot, tolerance = 1e-11)
    expect_equal(dl$measurement, cs$meas, tolerance = 1e-11)
    if (!is.null(cs$ep23)) {
      lps <- edge_probabilities(st, q_edges, log = TRUE)
      expect_equal(lps[2L], cs$ep23, tolerance = 1e-11)
    }
  }

  # LatentMultigraphBlockState
  lm_cases <- list(
    list(dir = FALSE, nest = FALSE, tot = 25.5780848648963),
    list(dir = FALSE, nest = TRUE,  tot = 29.654080288450206, dS01 = 0.44628710262842075),
    list(dir = TRUE,  nest = FALSE, tot = 30.501961069285276),
    list(dir = TRUE,  nest = TRUE,  tot = 34.89095637327525, dS01 = 3.1626361065453246)
  )
  for (cs in lm_cases) {
    g <- make_m10_ref_graph(cs$dir)
    st <- latent_multigraph_block_state(
      g, nested = cs$nest, blocks = if (cs$nest) bs_vals else b_vals
    )
    expect_equal(uncertain_description_length(st), cs$tot, tolerance = 1e-11)
    dSs <- rgraphtool:::cpp_uncertain_edge_dS(
      st$model_code, st$n_vertices, st$directed, st$self_loops,
      st$obs_edges, st$obs_v1, st$obs_v2, st$def1, st$def2, st$hparams,
      st$lp, st$lq, st$latent_edges, st$latent_weights, st$bs,
      st$degree_corrected, st$dist,
      matrix(c(1L, 2L, 1L, 4L), ncol = 2L, byrow = TRUE), c(1L, 1L), Inf
    )
    expect_true(is.infinite(dSs[2L]) && dSs[2L] > 0)
    if (!is.null(cs$dS01)) {
      expect_equal(dSs[1L], cs$dS01, tolerance = 1e-11)
    }
  }
})

test_that("latent_multigraph EM matches pinned graph-tool 3.9 weights", {
  g_u <- make_m10_ref_graph(FALSE)
  em_u <- latent_multigraph(g_u, epsilon = 1e-13)
  expect_s3_class(em_u, "rgraphtool_latent_multigraph")
  ku <- paste0(pmin(em_u$edges$from, em_u$edges$to), ":", pmax(em_u$edges$from, em_u$edges$to))
  wu <- setNames(em_u$edges$weight, ku)
  expect_equal(unname(wu["1:2"]), 1.1789473135222492, tolerance = 1e-10)
  expect_equal(unname(wu["3:4"]), 2.2673460919712456, tolerance = 1e-10)
  expect_equal(unname(wu["1:1"]), 0.169400098785355, tolerance = 1e-10)
  expect_equal(unname(wu["3:3"]), 0.9711263251744332, tolerance = 1e-10)
  expect_equal(unname(wu["7:7"]), 0.0, tolerance = 1e-12)

  g_d <- make_m10_ref_graph(TRUE)
  em_d <- latent_multigraph(g_d, epsilon = 1e-13)
  kd <- paste0(em_d$edges$from, ":", em_d$edges$to)
  wd <- setNames(em_d$edges$weight, kd)
  expect_equal(unname(wd["1:2"]), 1.07476892459312, tolerance = 1e-10)
  expect_equal(unname(wd["3:4"]), 1.7410491496063065, tolerance = 1e-10)
  expect_equal(unname(wd["1:1"]), 0.17217665316112074, tolerance = 1e-10)
  expect_equal(unname(wd["3:3"]), 0.5791056002597941, tolerance = 1e-10)
  expect_equal(unname(wd["7:7"]), 0.0, tolerance = 1e-12)
  expect_output(print(em_u), "rgraphtool_latent_multigraph")
})

test_that("incremental edge and block move deltas match full recomputation", {
  set.seed(101)
  b_vals <- c(1L, 1L, 1L, 2L, 2L, 2L, 3L)
  bs_vals <- list(b_vals, c(1L, 1L, 2L), c(1L, 1L))
  for (dir in c(FALSE, TRUE)) {
    g <- make_m10_ref_graph(dir)
    for (nest in c(FALSE, TRUE)) {
      bl <- if (nest) bs_vals else b_vals
      states <- list(
        measured_block_state(g, n = "n", x = "x", n_default = 2L, x_default = 0L, nested = nest, blocks = bl),
        mixed_measured_block_state(g, n = "n", x = "x", n_default = 2L, x_default = 0L, nested = nest, blocks = bl),
        uncertain_block_state(g, q = "q", q_default = 0.05, nested = nest, blocks = bl),
        latent_multigraph_block_state(g, nested = nest, blocks = bl)
      )
      for (st in states) {
        errs <- rgraphtool:::cpp_uncertain_check_moves(
          st$model_code, st$n_vertices, st$directed, st$self_loops,
          st$obs_edges, st$obs_v1, st$obs_v2, st$def1, st$def2, st$hparams,
          st$lp, st$lq, st$latent_edges, st$latent_weights, st$bs,
          st$degree_corrected, st$dist, 50L, Inf
        )
        expect_lt(max(errs), 1e-10)
      }
    }
  }
})

test_that("reconstruct_network recovers noisy measured SBM edges and error rates", {
  set.seed(2026)
  N <- 20L
  b_true <- rep(1:2, each = 10L)
  all_pairs <- t(utils::combn(N, 2L))
  same_blk <- b_true[all_pairs[, 1L]] == b_true[all_pairs[, 2L]]
  true_A <- as.integer(stats::runif(nrow(all_pairs)) < ifelse(same_blk, 0.65, 0.04))
  n_trials <- 3L
  p_fn <- 0.20
  q_fp <- 0.05
  x_obs <- stats::rbinom(nrow(all_pairs), n_trials, ifelse(true_A == 1L, 1 - p_fn, q_fp))
  keep <- x_obs > 0L
  g_obs <- graph_create(N, directed = FALSE)
  graph_add_edges(g_obs, all_pairs[keep, , drop = FALSE],
                  properties = list(n = rep(n_trials, sum(keep)), x = x_obs[keep]))

  st <- measured_block_state(g_obs, n = "n", x = "x", n_default = n_trials, x_default = 0L)
  expect_output(print(st), "rgraphtool_uncertain_block_state")

  rec <- reconstruct_network(st, nsamples = 120L, burnin = 100L, thin = 2L,
                             collect_partitions = TRUE)
  expect_s3_class(rec, "rgraphtool_reconstruction")
  expect_s3_class(rec$consensus_graph, "rgraphtool_graph")
  expect_length(rec$partitions, 120L)
  expect_output(print(rec), "rgraphtool_reconstruction")

  ek <- function(u, v) paste0(pmin(u, v), ":", pmax(u, v))
  idx <- match(ek(all_pairs[, 1L], all_pairs[, 2L]), ek(rec$edges$from, rec$edges$to))
  p_hat <- ifelse(is.na(idx), 0, rec$edges$prob[idx])
  ev <- evaluate_edge_prediction(p_hat, true_A)
  expect_gt(ev$roc_auc, 0.90)
  expect_lt(mean(rec$q_samples), 0.15)
  expect_lt(mean(rec$p_samples), 0.40)

  # MixedMeasured hyperparameter sampling smoke check
  st_mix <- mixed_measured_block_state(g_obs, n = "n", x = "x",
                                       n_default = n_trials, x_default = 0L, nested = FALSE)
  rec_mix <- reconstruct_network(st_mix, nsamples = 40L, burnin = 30L)
  expect_length(rec_mix$final_hparams, 4L)
  expect_true(all(rec_mix$final_hparams > 0))

  # Uncertain and LatentMultigraph reconstruction smoke checks
  g_ref <- make_m10_ref_graph(FALSE)
  rec_u <- reconstruct_network(uncertain_block_state(g_ref, q = "q", q_default = 0.05),
                               nsamples = 30L, burnin = 20L)
  expect_s3_class(rec_u, "rgraphtool_reconstruction")
  rec_lm <- reconstruct_network(latent_multigraph_block_state(g_ref),
                                nsamples = 30L, burnin = 20L)
  expect_true(all(rec_lm$edges$mean_weight > 0))
})

test_that("predict_edges and evaluate_edge_prediction work end-to-end", {
  set.seed(303)
  N <- 30L
  b_true <- rep(1:3, each = 10L)
  all_pairs <- t(utils::combn(N, 2L))
  same_blk <- b_true[all_pairs[, 1L]] == b_true[all_pairs[, 2L]]
  e_mask <- stats::runif(nrow(all_pairs)) < ifelse(same_blk, 0.75, 0.02)
  all_e <- all_pairs[e_mask, , drop = FALSE]
  hold_idx <- sample(nrow(all_e), 20L)
  g_train <- graph_create(N, directed = FALSE)
  graph_add_edges(g_train, all_e[-hold_idx, , drop = FALSE])

  held_out <- all_e[hold_idx, , drop = FALSE]
  non_e <- all_pairs[!e_mask, , drop = FALSE]
  neg_s <- non_e[sample(nrow(non_e), 20L), , drop = FALSE]
  cands <- rbind(held_out, neg_s)
  labels <- c(rep(1L, 20L), rep(0L, 20L))

  pred <- predict_edges(g_train, cands, nsamples = 40L, burnin = 30L, nested = TRUE)
  expect_s3_class(pred, "rgraphtool_edge_prediction")
  expect_equal(nrow(pred$predictions), 40L)
  expect_true(all(pred$predictions$prob >= 0 & pred$predictions$prob <= 1))
  expect_output(print(pred), "rgraphtool_edge_prediction")

  ev <- evaluate_edge_prediction(pred, labels)
  expect_s3_class(ev, "rgraphtool_edge_evaluation")
  expect_gt(ev$roc_auc, 0.80)
  expect_gt(ev$pr_auc, 0.70)
  expect_output(print(ev), "rgraphtool_edge_evaluation")

  # Flat block_state input
  fit_flat <- minimize_blockmodel_dl(g_train)
  pred_flat <- predict_edges(fit_flat, cands[1:6, , drop = FALSE], nsamples = 15L, burnin = 10L)
  expect_equal(pred_flat$model, "flat_sbm")
  expect_equal(nrow(pred_flat$predictions), 6L)

  # Uncertain state input
  g_ref <- make_m10_ref_graph(FALSE)
  st_u <- uncertain_block_state(g_ref, q = "q", q_default = 0.05)
  pred_u <- predict_edges(st_u, matrix(c(1L, 2L, 1L, 4L), ncol = 2L, byrow = TRUE),
                          nsamples = 20L, burnin = 10L)
  expect_equal(pred_u$model, "uncertain")

  # Exact evaluate_edge_prediction edge cases
  ev_perf <- evaluate_edge_prediction(c(0.9, 0.8, 0.2, 0.1), c(1, 1, 0, 0))
  expect_equal(ev_perf$roc_auc, 1.0)
  expect_equal(ev_perf$pr_auc, 1.0)
  expect_equal(ev_perf$f1, 1.0)

  ev_inv <- evaluate_edge_prediction(c(0.1, 0.2, 0.8, 0.9), c(1, 1, 0, 0))
  expect_equal(ev_inv$roc_auc, 0.0)

  ev_tie <- evaluate_edge_prediction(c(0.5, 0.5, 0.5, 0.5), c(1, 1, 0, 0))
  expect_equal(ev_tie$roc_auc, 0.5)
})

test_that("uncertain constructors validate inputs cleanly", {
  g <- make_m10_ref_graph(FALSE)
  expect_error(measured_block_state(g, n = rep(1, 7), x = rep(2, 7)), "0 <= x <= n")
  expect_error(measured_block_state(g, lp = 0.5), "non-positive")
  expect_error(measured_block_state(g, fn_params = list(alpha = -1, beta = 1)), "positive")
  expect_error(uncertain_block_state(g, q = rep(1.5, 7)), "\\[0, 1\\]")
  expect_error(uncertain_block_state(g, q_default = -0.1), "\\[0, 1\\]")
  expect_error(measured_error_posterior(uncertain_block_state(g, q = "q")), "measured_block_state")
  expect_error(evaluate_edge_prediction(c(0.2, 0.8), c(1, 1)), "at least one")
})
