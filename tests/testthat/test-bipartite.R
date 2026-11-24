test_that("bipartite_graph, is_bipartite, and as_incidence_matrix round-trip cleanly", {
  M <- matrix(
    c(1, 0, 2,
      0, 3, 0),
    nrow = 2L, byrow = TRUE,
    dimnames = list(plant = c("p1", "p2"), pollinator = c("a1", "a2", "a3"))
  )

  # Auto multigraph mode from integer counts
  g_multi <- bipartite_graph(M)
  expect_true(is_bipartite(g_multi))
  expect_equal(graph_num_vertices(g_multi), 5)
  expect_equal(graph_num_edges(g_multi), 6)
  expect_equal(graph_vertex_property(g_multi, "type"), c(FALSE, FALSE, TRUE, TRUE, TRUE))
  expect_equal(graph_vertex_property(g_multi, "name"), c("p1", "p2", "a1", "a2", "a3"))
  expect_equal(graph_graph_property(g_multi, "dim_labels"), c("plant", "pollinator"))

  M_back <- as_incidence_matrix(g_multi)
  expect_equal(M_back, M)

  # Weighted mode from continuous matrix
  Mw <- matrix(c(1.5, 0, -2.5, 0, 4.0, 0), nrow = 2L, byrow = TRUE)
  g_w <- bipartite_graph(Mw, weight_mode = "weight", dim_labels = c("u", "v"))
  expect_true(is_bipartite(g_w))
  expect_equal(graph_num_edges(g_w), 3)
  expect_equal(as_incidence_matrix(g_w, weight = "weight"), unname(Mw), ignore_attr = TRUE)

  # Reject non-bipartite graphs with within-type edges
  g_bad <- graph_create(3L, directed = FALSE)
  graph_add_edges(g_bad, rbind(c(1, 2), c(2, 3), c(1, 3)))
  expect_false(is_bipartite(g_bad))
  expect_error(bipartite_graph(g_bad, type = c(FALSE, FALSE, TRUE)), "not bipartite")

  # Round-trip through .gt and igraph
  f <- tempfile(fileext = ".gt")
  on.exit(unlink(f), add = TRUE)
  write_gt(g_multi, f)
  g_gt <- read_gt(f)
  expect_true(is_bipartite(g_gt))
  expect_equal(as_incidence_matrix(g_gt), M)

  if (requireNamespace("igraph", quietly = TRUE)) {
    ig <- as_igraph(g_multi)
    expect_true(igraph::is_bipartite(ig))
    g_ig <- bipartite_graph(ig)
    expect_equal(as_incidence_matrix(g_ig), M)
  }
})

test_that("bipartite_description_length matches pinned graph-tool 3.9 pclabel to < 1e-12", {
  edges <- rbind(c(1, 5), c(1, 5), c(1, 6), c(2, 5), c(2, 7), c(3, 6), c(3, 7))
  blocks <- c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L)
  type <- c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L)

  # Undirected graph (7 edges, 2 isolated vertices)
  gu <- graph_create(8L, directed = FALSE)
  graph_add_edges(gu, edges)

  expect_equal(
    bipartite_description_length(bipartite_block_state(
      gu, type = type, blocks = blocks, degree_corrected = TRUE,
      degree_prior = "distributed", edge_prior = "unipartite"
    )),
    26.374574699861881, tolerance = 1e-12
  )
  expect_equal(
    bipartite_description_length(bipartite_block_state(
      gu, type = type, blocks = blocks, degree_corrected = TRUE,
      degree_prior = "uniform", edge_prior = "unipartite"
    )),
    25.799210554958321, tolerance = 1e-12
  )
  expect_equal(
    bipartite_description_length(bipartite_block_state(
      gu, type = type, blocks = blocks, degree_corrected = FALSE,
      edge_prior = "unipartite"
    )),
    25.117357381017133, tolerance = 1e-12
  )

  # Exact Yen & Larremore (2020) bipartite edge_prior difference:
  # B0 = 2, B1 = 2, B = 4, E = 7
  s_uni <- bipartite_description_length(bipartite_block_state(
    gu, type = type, blocks = blocks, degree_corrected = TRUE, edge_prior = "unipartite"
  ))
  s_bip <- bipartite_description_length(bipartite_block_state(
    gu, type = type, blocks = blocks, degree_corrected = TRUE, edge_prior = "bipartite"
  ))
  expected_diff <- lchoose(2 * 2 + 7 - 1, 7) - lchoose(4 * 5 / 2 + 7 - 1, 7)
  expect_equal(s_bip - s_uni, expected_diff, tolerance = 1e-13)

  # Directed bipartite graph (8 edges including reverse 6 -> 2)
  gd <- graph_create(8L, directed = TRUE)
  graph_add_edges(gd, rbind(edges, c(6L, 2L)))

  expect_equal(
    bipartite_description_length(bipartite_block_state(
      gd, type = type, blocks = blocks, degree_corrected = TRUE,
      degree_prior = "distributed", edge_prior = "unipartite"
    )),
    31.518799077226220, tolerance = 1e-12
  )
  expect_equal(
    bipartite_description_length(bipartite_block_state(
      gd, type = type, blocks = blocks, degree_corrected = TRUE,
      degree_prior = "uniform", edge_prior = "unipartite"
    )),
    30.943434932322660, tolerance = 1e-12
  )
  expect_equal(
    bipartite_description_length(bipartite_block_state(
      gd, type = type, blocks = blocks, degree_corrected = FALSE,
      edge_prior = "unipartite"
    )),
    30.261581758381475, tolerance = 1e-12
  )
})

test_that("incremental move_delta matches full entropy recomputation to < 1e-12", {
  set.seed(123)
  edges <- rbind(
    c(1, 6), c(1, 6), c(1, 7), c(2, 6), c(2, 8),
    c(3, 7), c(3, 9), c(4, 8), c(4, 9), c(5, 10)
  )
  type <- c(rep(0L, 5L), rep(1L, 5L))
  init_b <- c(1L, 1L, 2L, 2L, 3L, 4L, 4L, 5L, 5L, 6L)

  for (directed in c(FALSE, TRUE)) {
    g <- graph_create(10L, directed = directed)
    e_use <- if (directed) rbind(edges, c(7L, 2L), c(9L, 4L)) else edges
    w <- seq(0.5, by = 0.7, length.out = nrow(e_use))
    graph_add_edges(g, e_use, properties = list(w = w))

    for (dc in c(TRUE, FALSE)) {
      for (dist in (if (dc) c(TRUE, FALSE) else TRUE)) {
        for (bp_prior in c(TRUE, FALSE)) {
          diffs <- rgraphtool:::cpp_bipartite_check_moves(
            g$pointer, type, init_b, dc, dist, bp_prior,
            matrix(0, 0, 0), character(0), list(), 100L, 1000
          )
          expect_lt(max(abs(diffs)), 1e-12)
        }
      }
    }

    # Weighted incremental checks (real-exponential & real-normal)
    W <- matrix(graph_edge_property(g, "w"), ncol = 1L)
    diffs_exp <- rgraphtool:::cpp_bipartite_check_moves(
      g$pointer, type, init_b, TRUE, TRUE, TRUE,
      W, "real-exponential", list(c(NaN, NaN)), 100L, 1000
    )
    expect_lt(max(abs(diffs_exp)), 1e-12)

    diffs_norm <- rgraphtool:::cpp_bipartite_check_moves(
      g$pointer, type, init_b, TRUE, TRUE, TRUE,
      W, "real-normal", list(c(mean(w), 1, mean((w - mean(w))^2), 3)), 100L, 1000
    )
    expect_lt(max(abs(diffs_norm)), 1e-12)
  }
})

test_that("minimize_bipartite_blockmodel_dl recovers planted binary, weighted, and nested blocks", {
  set.seed(42)

  # 1. Binary planted 2 x 2 bipartite SBM
  M <- matrix(0L, 24L, 24L)
  M[1:12, 1:12] <- stats::rbinom(144, 1, 0.85)
  M[13:24, 13:24] <- stats::rbinom(144, 1, 0.85)
  M[1:12, 13:24] <- stats::rbinom(144, 1, 0.04)
  M[13:24, 1:12] <- stats::rbinom(144, 1, 0.04)

  fit <- minimize_bipartite_blockmodel_dl(M, dim_labels = c("row", "col"), nruns = 2L)
  expect_s3_class(fit, "rgraphtool_bipartite_block_state")
  expect_equal(unname(fit$B), c(2L, 2L))
  mem <- bipartite_membership(fit)
  expect_equal(names(mem), c("row", "col"))
  expect_equal(partition_overlap(mem$row, rep(1:2, each = 12L)), 1.0)
  expect_equal(partition_overlap(mem$col, rep(1:2, each = 12L)), 1.0)

  # Respect B_row_range / B_col_range
  fit_c <- minimize_bipartite_blockmodel_dl(M, B_row_range = c(3L, 4L), B_col_range = c(2L, 2L))
  expect_true(fit_c$B[1L] >= 3L && fit_c$B[1L] <= 4L)
  expect_equal(unname(fit_c$B[2L]), 2L)

  # 2. Weighted bipartite SBM where topology is uniform complete, only weights carry blocks
  Mw <- matrix(0, 16L, 16L)
  Mw[1:8, 1:8] <- stats::rnorm(64, mean = 10, sd = 0.5)
  Mw[9:16, 9:16] <- stats::rnorm(64, mean = 10, sd = 0.5)
  Mw[1:8, 9:16] <- stats::rnorm(64, mean = -10, sd = 0.5)
  Mw[9:16, 1:8] <- stats::rnorm(64, mean = -10, sd = 0.5)
  fit_w <- minimize_bipartite_blockmodel_dl(
    Mw, weight_mode = "weight", weights = "weight", types = "real-normal", nruns = 2L
  )
  expect_equal(unname(fit_w$B), c(2L, 2L))
  mw_mat <- bipartite_block_matrix(fit_w, normalize = "mean_weight")
  expect_equal(dim(mw_mat), c(2L, 2L))

  # 3. Nested bipartite SBM & model comparison
  nfit <- minimize_nested_bipartite_blockmodel_dl(M, dim_labels = c("row", "col"), nruns = 2L)
  expect_s3_class(nfit, "rgraphtool_nested_bipartite_block_state")
  expect_true(is.finite(nested_bipartite_description_length(nfit)))
  expect_equal(unname(nfit$B), c(2L, 2L))

  cmp <- compare_bipartite_models(M, nruns = 1L)
  expect_s3_class(cmp, "rgraphtool_bipartite_model_comparison")
  expect_equal(nrow(cmp$comparison), 4L)
})

test_that("bipartite MCMC posterior sampling and plotting work cleanly", {
  set.seed(7)
  M <- matrix(0L, 12L, 12L)
  M[1:6, 1:6] <- stats::rbinom(36, 1, 0.8)
  M[7:12, 7:12] <- stats::rbinom(36, 1, 0.8)
  M[1:6, 7:12] <- stats::rbinom(36, 1, 0.1)
  M[7:12, 1:6] <- stats::rbinom(36, 1, 0.1)

  fit <- minimize_bipartite_blockmodel_dl(M, dim_labels = c("plant", "pollinator"))
  sw <- mcmc_sweep(fit, beta = 1, niter = 10L)
  expect_s3_class(sw$state, "rgraphtool_bipartite_block_state")
  expect_true(is.finite(sw$entropy))

  post <- sample_posterior(fit, nsamples = 15L, burnin = 10L, thin = 2L)
  expect_length(post$samples, 15L)
  # Every sampled partition must keep type-0 and type-1 blocks disjoint
  for (s in post$samples) {
    expect_length(intersect(unique(s[fit$type == 0L]), unique(s[fit$type == 1L])), 0L)
  }

  # Plotting all 3 modes to a temporary PDF
  pdf_file <- tempfile(fileext = ".pdf")
  on.exit(unlink(pdf_file), add = TRUE)
  grDevices::pdf(pdf_file)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot_bipartite(fit, type = "matrix"))
  expect_silent(plot_bipartite(fit, type = "meso"))
  expect_silent(plot_bipartite(fit, type = "network"))
  expect_true(file.exists(pdf_file))
})

test_that("bipartite edge cases (empty graph, E = 0) are handled cleanly", {
  g0 <- graph_create(0L, directed = FALSE)
  st0 <- bipartite_block_state(g0, type = logical(0))
  expect_equal(bipartite_description_length(st0), 0)

  M0 <- matrix(0L, 4L, 5L)
  fit0 <- minimize_bipartite_blockmodel_dl(M0)
  expect_equal(graph_num_edges(fit0$graph), 0)
  expect_true(is.finite(bipartite_description_length(fit0)))
})
