test_that("flat and nested label alignment and overlap match graph-tool 3.9", {
  x <- c(1L, 1L, 1L, 2L, 2L, 2L)
  y <- c(3L, 3L, 2L, 2L, 2L, 1L)
  expect_equal(partition_overlap(x, y, norm = TRUE), 4 / 6, tolerance = 1e-12)
  expect_equal(partition_overlap(x, y, norm = FALSE), 4, tolerance = 1e-12)
  expect_equal(align_partition_labels(y, x), c(1L, 1L, 2L, 2L, 2L, 3L))

  # Permuting labels preserves overlap = 1 and recovers original partition
  perm <- c(4L, 4L, 4L, 2L, 2L, 2L)
  expect_equal(partition_overlap(perm, x), 1, tolerance = 1e-12)
  expect_equal(align_partition_labels(perm, x), x)

  # Nested alignment and overlap (graph-tool 3.9 reference: 0.625 and 8)
  nx <- list(c(1L, 1L, 1L, 2L, 2L, 2L), c(1L, 2L), c(1L, 1L))
  ny <- list(c(2L, 2L, 1L, 1L, 1L, 3L), c(2L, 1L, 2L), c(1L, 1L))
  expect_equal(nested_partition_overlap(nx, ny, norm = TRUE), 0.625, tolerance = 1e-12)
  expect_equal(nested_partition_overlap(nx, ny, norm = FALSE), 8, tolerance = 1e-12)
  expect_equal(nested_partition_overlap(nx, nx, norm = TRUE), 1, tolerance = 1e-12)
})

test_that("PartitionModeState and ModeClusterState scores match graph-tool 3.9 to 1e-10", {
  bs_fixed <- list(
    c(1L, 1L, 1L, 2L, 2L, 2L),
    c(1L, 1L, 1L, 2L, 2L, 2L),
    c(1L, 1L, 2L, 2L, 2L, 3L),
    c(1L, 1L, 1L, 2L, 2L, 3L)
  )
  pm <- partition_mode_state(bs_fixed, relabel = FALSE, nested = FALSE, converge = FALSE)
  expect_equal(pm$entropy, 21.21811450618926, tolerance = 1e-10)
  expect_equal(mode_entropy(pm, mle = TRUE), 1.2554823251787535, tolerance = 1e-10)
  expect_equal(mode_entropy(pm, mle = FALSE), 5.145188924861863, tolerance = 1e-10)
  expect_equal(mode_uncertainty(pm, mle = TRUE), 0.125, tolerance = 1e-10)
  expect_equal(mode_uncertainty(pm, mle = FALSE), 5 / 14, tolerance = 1e-10)
  expect_equal(rowSums(mode_marginals(pm)), rep(4L, 6L))
  expect_equal(rowSums(mode_marginals(pm, probability = TRUE)), rep(1, 6L), tolerance = 1e-12)
  expect_equal(mode_consensus(pm), c(1L, 1L, 1L, 2L, 2L, 2L))

  mc <- mode_cluster_state(bs_fixed, nested = FALSE, b = c(1L, 1L, 2L, 2L),
                           relabel = FALSE, niter = 0L, converge = FALSE)
  expect_equal(mc$entropy, 25.895562785409098, tolerance = 1e-10)
  expect_equal(mc$K, 2L)
  expect_equal(mc$Be, 2, tolerance = 1e-10)
  expect_equal(mc$weights, c(0.5, 0.5))

  # Nested PartitionModeState (graph-tool 3.9 reference)
  nx <- list(c(1L, 1L, 1L, 2L, 2L, 2L), c(1L, 2L), c(1L, 1L))
  ny <- list(c(2L, 2L, 1L, 1L, 1L, 3L), c(2L, 1L, 2L), c(1L, 1L))
  pmn <- partition_mode_state(list(nx, ny), relabel = TRUE, nested = TRUE, converge = TRUE)
  expect_equal(pmn$entropy, 18.610741872840492, tolerance = 1e-10)
  expect_equal(mode_uncertainty(pmn, mle = TRUE), 1 / 6, tolerance = 1e-10)
  expect_equal(mode_entropy(pmn, mle = TRUE), log(4), tolerance = 1e-10)
  expect_equal(dim(mode_marginals(pmn, level = 0L)), c(6L, 3L))
})

test_that("mode_cluster_state separates distinct synthetic posterior modes", {
  set.seed(42)
  n <- 40L
  # Mode A: two halves (1..20, 21..40)
  # Mode B: alternating blocks (odd vs even)
  base_a <- rep(1:2, each = 20L)
  base_b <- rep(1:2, times = 20L)
  make_noisy <- function(base, flips = 1L) {
    b <- base
    idx <- sample.int(n, flips)
    b[idx] <- 3L - b[idx]
    # randomly permute labels to verify label-invariance
    perm <- sample.int(2L)
    perm[b]
  }
  parts <- c(
    replicate(30L, make_noisy(base_a, 1L), simplify = FALSE),
    replicate(20L, make_noisy(base_b, 1L), simplify = FALSE)
  )
  parts <- parts[sample.int(length(parts))]

  mc <- mode_cluster_state(parts, nested = FALSE, niter = 20L)
  expect_equal(mc$K, 2L)
  expect_equal(sort(vapply(mc$modes, function(m) m$M, integer(1))), c(20L, 30L))
  expect_gt(partition_overlap(mode_consensus(mc$modes[[1L]]), base_a), 0.95)
  expect_gt(partition_overlap(mode_consensus(mc$modes[[2L]]), base_b), 0.95)
  expect_lt(mode_uncertainty(mc$modes[[1L]]), 0.1)
  expect_output(print(mc), "mode cluster state")
  expect_output(print(mc$modes[[1L]]), "partition mode")
})

test_that("posterior sampling integrates cleanly with partition_mode_state and mode_cluster_state", {
  set.seed(55)
  g <- graph_create(36, directed = FALSE)
  edges <- rbind(
    do.call(rbind, lapply(list(1:12, 13:24, 25:36), function(v) t(combn(v, 2)))),
    c(1, 13), c(13, 25), c(25, 1)
  )
  graph_add_edges(g, edges)
  fit <- minimize_blockmodel_dl(g)
  post <- sample_posterior(fit, nsamples = 60L, thin = 2L, burnin = 20L)
  pm <- partition_mode_state(post)
  mc <- mode_cluster_state(post)
  expect_equal(pm$M, 60L)
  expect_gte(mc$K, 1L)
  expect_equal(sum(vapply(mc$modes, function(m) m$M, integer(1))), 60L)
  expect_gt(partition_overlap(mode_consensus(pm), rep(1:3, each = 12L)), 0.9)

  # Nested posterior sampling -> nested mode clustering
  nfit <- minimize_nested_blockmodel_dl(g)
  npost <- sample_posterior(nfit, nsamples = 40L, thin = 2L, burnin = 15L)
  npm <- partition_mode_state(npost)
  nmc <- mode_cluster_state(npost)
  expect_true(npm$nested)
  expect_equal(npm$M, 40L)
  expect_true(nmc$nested)
  expect_equal(nrow(mode_summary(nmc)), nmc$K)
})
