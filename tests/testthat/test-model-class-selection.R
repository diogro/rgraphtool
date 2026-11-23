test_that("mode_lprob, cluster_lprob, and nested/cluster mode_entropy match graph-tool 3.9", {
  # Flat partition mode reference from graph-tool 3.9
  bs <- list(
    c(1L, 1L, 1L, 2L, 2L, 2L),
    c(1L, 1L, 1L, 2L, 2L, 2L),
    c(1L, 1L, 2L, 2L, 2L, 3L),
    c(1L, 1L, 1L, 2L, 2L, 3L)
  )
  pm <- partition_mode_state(bs, relabel = FALSE, converge = FALSE)
  expect_equal(mode_lprob(pm, c(1L, 1L, 1L, 2L, 2L, 2L), mle = TRUE),
               -0.9808292530117261, tolerance = 1e-12)
  expect_equal(mode_lprob(pm, c(1L, 1L, 1L, 2L, 2L, 2L), mle = FALSE),
               -2.752802594807478, tolerance = 1e-12)
  expect_equal(mode_lprob(pm, c(1L, 1L, 2L, 2L, 2L, 3L), mle = TRUE),
               -2.0794415416798357, tolerance = 1e-12)
  expect_equal(mode_lprob(pm, c(1L, 1L, 2L, 2L, 2L, 3L), mle = FALSE),
               -3.4459497753674233, tolerance = 1e-12)
  # Misaligned partition without align = TRUE vs with align = TRUE
  expect_equal(mode_lprob(pm, c(2L, 2L, 2L, 1L, 1L, 1L), mle = TRUE, align = FALSE), -Inf)
  expect_equal(mode_lprob(pm, c(2L, 2L, 2L, 1L, 1L, 1L), mle = FALSE, align = FALSE),
               -10.982313713771935, tolerance = 1e-12)
  expect_equal(mode_lprob(pm, c(2L, 2L, 2L, 1L, 1L, 1L), mle = TRUE, align = TRUE),
               -0.9808292530117261, tolerance = 1e-12)
  # Unseen label returns -Inf for both MLE and posterior mean
  expect_equal(mode_lprob(pm, c(1L, 1L, 1L, 2L, 2L, 6L), mle = TRUE), -Inf)
  expect_equal(mode_lprob(pm, c(1L, 1L, 1L, 2L, 2L, 6L), mle = FALSE), -Inf)

  # ModeClusterState reference from graph-tool 3.9
  bs2 <- list(
    c(1L, 1L, 1L, 2L, 2L, 2L),
    c(1L, 1L, 1L, 2L, 2L, 2L),
    c(1L, 1L, 1L, 2L, 2L, 2L),
    c(1L, 2L, 1L, 2L, 1L, 2L),
    c(1L, 2L, 1L, 2L, 1L, 2L)
  )
  mc <- mode_cluster_state(bs2, b = c(1L, 1L, 1L, 2L, 2L), niter = 0L,
                           relabel = FALSE, converge = FALSE)
  expect_equal(mode_entropy(mc, mle = TRUE), 0.6730116670092565, tolerance = 1e-12)
  expect_equal(mode_entropy(mc, mle = FALSE), 3.824064738831873, tolerance = 1e-12)
  expect_equal(cluster_lprob(mc, c(1L, 1L, 1L, 2L, 2L, 2L), mode = 1L, mle = TRUE),
               log(3 / 5), tolerance = 1e-12)
  expect_equal(cluster_lprob(mc, c(1L, 1L, 1L, 2L, 2L, 2L), mode = 1L, mle = FALSE),
               -1.84968693, tolerance = 1e-6)
  expect_equal(cluster_lprob(mc, c(1L, 2L, 1L, 2L, 1L, 2L), mode = 2L, mle = TRUE),
               log(2 / 5), tolerance = 1e-12)
  expect_equal(cluster_lprob(mc, c(1L, 2L, 1L, 2L, 1L, 2L), mode = 2L, mle = FALSE),
               -2.64238317, tolerance = 1e-6)
  # Full mixture log-prob (mode = NULL) equals log(3/5) when only mode 1 supports partition
  expect_equal(cluster_lprob(mc, c(1L, 1L, 1L, 2L, 2L, 2L), mle = TRUE),
               log(3 / 5), tolerance = 1e-12)

  # Two-mode cluster with within-mode uncertainty
  bs3 <- list(
    c(1L, 1L, 1L, 2L, 2L, 2L),
    c(1L, 1L, 2L, 2L, 2L, 2L),
    c(1L, 1L, 1L, 2L, 2L, 2L),
    c(1L, 2L, 1L, 2L, 1L, 2L),
    c(1L, 2L, 1L, 2L, 1L, 1L)
  )
  mc3 <- mode_cluster_state(bs3, b = c(1L, 1L, 1L, 2L, 2L), niter = 0L,
                            relabel = FALSE, converge = FALSE)
  expect_equal(mode_entropy(mc3, mle = TRUE), 1.3321790402101223, tolerance = 1e-12)

  # Nested PartitionModeState posterior_entropy & posterior_lprob reference from graph-tool 3.9
  nbs <- list(
    list(c(1L, 1L, 1L, 2L, 2L, 2L), c(1L, 1L)),
    list(c(1L, 1L, 1L, 2L, 2L, 2L), c(1L, 2L), c(1L, 1L)),
    list(c(1L, 1L, 2L, 2L, 2L, 2L), c(1L, 1L))
  )
  npm <- partition_mode_state(nbs, nested = TRUE, relabel = FALSE, converge = FALSE)
  expect_equal(mode_entropy(npm, mle = TRUE), 1.273028336589626, tolerance = 1e-12)
  expect_equal(mode_entropy(npm, mle = FALSE), 4.34843787524764, tolerance = 1e-12)
  expect_equal(mode_lprob(npm, nbs[[1L]], mle = TRUE), -0.810930216216, tolerance = 1e-10)
  expect_equal(mode_lprob(npm, nbs[[1L]], mle = FALSE), -2.360512555417, tolerance = 1e-10)
  expect_equal(mode_lprob(npm, nbs[[2L]], mle = TRUE), -1.504077396776, tolerance = 1e-10)
  expect_equal(mode_lprob(npm, nbs[[2L]], mle = FALSE), -2.765977663525, tolerance = 1e-10)
})

test_that("model_evidence works across flat, nested, overlap, and nested overlap states", {
  # Two clean 6-node cliques connected by 1 bridge edge (sharp posterior)
  g <- graph_create(n = 12L, directed = FALSE)
  ed <- list()
  for (i in 1:5) for (j in (i + 1):6) ed[[length(ed) + 1L]] <- c(i, j)
  for (i in 7:11) for (j in (i + 1):12) ed[[length(ed) + 1L]] <- c(i, j)
  ed[[length(ed) + 1L]] <- c(6L, 7L)
  graph_add_edges(g, do.call(rbind, ed))

  set.seed(42)
  st_flat <- minimize_blockmodel_dl(g, degree_corrected = TRUE)
  ev_flat <- model_evidence(st_flat, nsamples = 40L, thin = 2L, burnin = 10L)
  expect_s3_class(ev_flat, "graphtool_model_evidence")
  expect_true(is.finite(ev_flat$log_evidence))
  expect_true(is.finite(ev_flat$importance_log_evidence))
  expect_equal(ev_flat$neg_log_evidence, -ev_flat$log_evidence)
  expect_gte(ev_flat$posterior_entropy, 0)
  # Variational and importance log-evidence should agree closely on sharp posterior
  expect_lt(abs(ev_flat$variational_log_evidence - ev_flat$importance_log_evidence), 2.5)
  expect_output(print(ev_flat), "model evidence")

  # Pre-collected samples/entropies interface
  ev_pre <- model_evidence(st_flat, samples = ev_flat$samples,
                           entropies = ev_flat$entropies, method = "importance")
  expect_equal(ev_pre$log_evidence, ev_flat$importance_log_evidence, tolerance = 1e-10)

  # Nested non-overlap
  set.seed(42)
  st_nested <- minimize_nested_blockmodel_dl(g, degree_corrected = TRUE)
  ev_nested <- model_evidence(st_nested, nsamples = 30L, thin = 2L, burnin = 10L)
  expect_true(is.finite(ev_nested$log_evidence))
  expect_true(is.finite(ev_nested$importance_log_evidence))

  # Flat overlap
  set.seed(42)
  st_ov <- minimize_overlap_blockmodel_dl(g, degree_corrected = TRUE)
  ev_ov <- model_evidence(st_ov, nsamples = 25L, thin = 2L, burnin = 10L)
  expect_true(is.finite(ev_ov$log_evidence))
  expect_true(is.finite(ev_ov$importance_log_evidence))

  # Nested overlap
  set.seed(42)
  st_nov <- minimize_nested_overlap_blockmodel_dl(g, degree_corrected = TRUE)
  ev_nov <- model_evidence(st_nov, nsamples = 20L, thin = 1L, burnin = 5L)
  expect_true(is.finite(ev_nov$log_evidence))
  expect_true(is.finite(ev_nov$importance_log_evidence))

  # Empty graph edge case for overlap
  g0 <- graph_create(n = 3L, directed = FALSE)
  ev0 <- model_evidence(overlap_block_state(g0), nsamples = 10L, burnin = 0L)
  expect_equal(ev0$log_evidence, 0)
})

test_that("compare_model_classes ranks all 8 models and computes Bayes factors", {
  # Two 6-node cliques with bridge
  g <- graph_create(n = 12L, directed = FALSE)
  ed <- list()
  for (i in 1:5) for (j in (i + 1):6) ed[[length(ed) + 1L]] <- c(i, j)
  for (i in 7:11) for (j in (i + 1):12) ed[[length(ed) + 1L]] <- c(i, j)
  ed[[length(ed) + 1L]] <- c(6L, 7L)
  graph_add_edges(g, do.call(rbind, ed))

  set.seed(42)
  cmp_map <- compare_model_classes(g, method = "map", nruns = 1L)
  expect_s3_class(cmp_map, "graphtool_model_class_comparison")
  expect_equal(nrow(cmp_map$summary), 8L)
  expect_equal(cmp_map$summary$delta[1L], 0)
  expect_true(all(diff(cmp_map$summary$neg_log_evidence) >= -1e-12))
  expect_equal(sum(cmp_map$summary$posterior_weight), 1, tolerance = 1e-12)
  expect_equal(cmp_map$best_model, cmp_map$summary$model[1L])
  # Log Bayes factor matrix is antisymmetric with zero diagonal
  lbf <- cmp_map$log_bayes_factors
  expect_equal(diag(lbf), rep(0, 8), ignore_attr = TRUE)
  expect_equal(lbf, -t(lbf), tolerance = 1e-12)
  expect_output(print(cmp_map), "model class comparison")

  # Posterior evidence comparison on subset of models
  set.seed(42)
  cmp_post <- compare_model_classes(
    g, models = c("flat_ndc", "flat_dc", "nested_ndc", "nested_dc"),
    method = "posterior", nsamples = 30L, thin = 2L, burnin = 10L
  )
  expect_equal(nrow(cmp_post$summary), 4L)
  expect_true(all(is.finite(cmp_post$summary$neg_log_evidence)))
  expect_true(all(is.finite(cmp_post$summary$posterior_entropy)))
  expect_equal(sum(cmp_post$summary$posterior_weight), 1, tolerance = 1e-12)
  expect_output(print(cmp_post), "posterior/variational")
})

test_that("compare_model_classes recovers true generative model class on synthetic graphs", {
  # 1. Strong degree heterogeneity (two star-hub communities -> favors degree-corrected!)
  g_hub <- graph_create(n = 30L, directed = FALSE)
  ed <- list()
  # Community 1: vertex 1 is a hub connected to 2..15, plus a sparse ring
  for (v in 2:15) ed[[length(ed) + 1L]] <- c(1L, v)
  for (v in 2:14) ed[[length(ed) + 1L]] <- c(v, v + 1L)
  # Community 2: vertex 16 is a hub connected to 17..30, plus a sparse ring
  for (v in 17:30) ed[[length(ed) + 1L]] <- c(16L, v)
  for (v in 17:29) ed[[length(ed) + 1L]] <- c(v, v + 1L)
  ed[[length(ed) + 1L]] <- c(1L, 16L)
  graph_add_edges(g_hub, do.call(rbind, ed))

  set.seed(42)
  cmp_hub <- compare_model_classes(g_hub, models = c("flat_ndc", "flat_dc"), method = "map")
  expect_equal(cmp_hub$best_model, "flat_dc")

  # 2. Pairwise overlapping 4-community network -> favors overlap over flat non-overlap!
  set.seed(123)
  cores <- list(1:8, 9:16, 17:24, 25:32)
  pairs <- utils::combn(4L, 2L, simplify = FALSE)
  groups <- cores
  idx <- 33L
  for (k in seq_along(pairs)) {
    br <- idx:(idx + 2L)
    idx <- idx + 3L
    groups[[pairs[[k]][1L]]] <- c(groups[[pairs[[k]][1L]]], br)
    groups[[pairs[[k]][2L]]] <- c(groups[[pairs[[k]][2L]]], br)
  }
  ed_ov <- list()
  for (cl in groups) {
    for (a in seq_along(cl)) {
      for (b in seq_along(cl)) {
        if (a < b && stats::runif(1) < 0.65) {
          ed_ov[[length(ed_ov) + 1L]] <- sort(c(cl[a], cl[b]))
        }
      }
    }
  }
  g_ov <- graph_create(n = 50L, directed = FALSE)
  graph_add_edges(g_ov, do.call(rbind, ed_ov))

  set.seed(42)
  cmp_ov <- compare_model_classes(
    g_ov, models = c("flat_ndc", "flat_dc", "overlap_flat_ndc", "overlap_flat_dc"),
    method = "map"
  )
  expect_true(grepl("^overlap_", cmp_ov$best_model))

  # Validation errors
  expect_error(compare_model_classes(g_hub, models = "not_a_model"), "unknown model")
  expect_error(compare_model_classes(g_hub, nruns = 0L), "nruns")
  expect_error(model_evidence(list()), "state must be")
})
