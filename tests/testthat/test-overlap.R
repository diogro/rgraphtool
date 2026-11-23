test_that("flat and nested OverlapBlockState entropies match graph-tool 3.9 on complex undirected graph", {
  # Undirected graph with parallel edges, self-loops, and isolated vertex (vertex 8)
  g <- graph_create(n = 8L, directed = FALSE)
  edges <- matrix(c(
    1L, 2L,  1L, 2L,  1L, 3L,  2L, 3L,  2L, 2L,  3L, 4L,  3L, 7L,
    4L, 5L,  4L, 6L,  4L, 6L,  5L, 6L,  6L, 6L,  6L, 6L,  7L, 7L
  ), ncol = 2L, byrow = TRUE)
  graph_add_edges(g, edges)
  he <- matrix(c(
    1L, 1L,  1L, 2L,  1L, 1L,  1L, 2L,  1L, 1L,  2L, 3L,  1L, 3L,
    3L, 3L,  3L, 3L,  3L, 3L,  3L, 2L,  2L, 2L,  2L, 3L,  3L, 3L
  ), ncol = 2L, byrow = TRUE)

  # Reference values from graph-tool 3.9 OverlapBlockState / NestedBlockState
  st_ndc <- overlap_block_state(g, half_edges = he, degree_corrected = FALSE)
  st_dc_dist <- overlap_block_state(g, half_edges = he, degree_corrected = TRUE,
                                    degree_prior = "distributed")
  st_dc_unif <- overlap_block_state(g, half_edges = he, degree_corrected = TRUE,
                                    degree_prior = "uniform")
  expect_equal(overlap_description_length(st_ndc), 44.3980815475771, tolerance = 1e-10)
  expect_equal(overlap_description_length(st_dc_dist), 46.98647957149643, tolerance = 1e-10)
  expect_equal(overlap_description_length(st_dc_unif), 48.796588179392685, tolerance = 1e-10)

  # Nested overlap (3 levels: level 0 = he, level 1 = [1, 2, 1], level 2 = [1, 1])
  nst_ndc <- nested_overlap_block_state(
    g, half_edges = he, upper_levels = list(c(1L, 2L, 1L), c(1L, 1L)),
    degree_corrected = FALSE
  )
  nst_dc <- nested_overlap_block_state(
    g, half_edges = he, upper_levels = list(c(1L, 2L, 1L), c(1L, 1L)),
    degree_corrected = TRUE
  )
  expect_equal(nested_overlap_description_length(nst_ndc), 49.00634292613486, tolerance = 1e-10)
  expect_equal(nested_overlap_description_length(nst_dc), 51.5947409500542, tolerance = 1e-10)
  expect_equal(sum(nested_overlap_level_entropy(nst_dc)),
               nested_overlap_description_length(nst_dc), tolerance = 1e-12)
})

test_that("flat and nested OverlapBlockState entropies match graph-tool 3.9 on complex directed graph", {
  # Directed graph with parallel edges, self-loops, reciprocal edges, and isolated vertex (vertex 7)
  g <- graph_create(n = 7L, directed = TRUE)
  edges <- matrix(c(
    1L, 2L,  1L, 2L,  1L, 3L,  2L, 1L,  2L, 2L,  3L, 4L,  3L, 5L,
    4L, 5L,  4L, 6L,  5L, 4L,  5L, 6L,  6L, 6L,  6L, 6L
  ), ncol = 2L, byrow = TRUE)
  graph_add_edges(g, edges)
  # Oriented as [source_block, target_block] matching graph-tool's min/max storage
  he <- matrix(c(
    1L, 2L,  1L, 2L,  1L, 1L,  1L, 2L,  2L, 2L,  1L, 3L,  1L, 2L,
    3L, 2L,  3L, 1L,  3L, 2L,  2L, 1L,  1L, 1L,  1L, 2L
  ), ncol = 2L, byrow = TRUE)

  st_ndc <- overlap_block_state(g, half_edges = he, degree_corrected = FALSE)
  st_dc_dist <- overlap_block_state(g, half_edges = he, degree_corrected = TRUE,
                                    degree_prior = "distributed")
  st_dc_unif <- overlap_block_state(g, half_edges = he, degree_corrected = TRUE,
                                    degree_prior = "uniform")
  expect_equal(overlap_description_length(st_ndc), 52.428622836954666, tolerance = 1e-10)
  expect_equal(overlap_description_length(st_dc_dist), 51.27808252947775, tolerance = 1e-10)
  expect_equal(overlap_description_length(st_dc_unif), 56.68972858133279, tolerance = 1e-10)

  nst_ndc <- nested_overlap_block_state(
    g, half_edges = he, upper_levels = list(c(1L, 2L, 1L), c(1L, 1L)),
    degree_corrected = FALSE
  )
  nst_dc <- nested_overlap_block_state(
    g, half_edges = he, upper_levels = list(c(1L, 2L, 1L), c(1L, 1L)),
    degree_corrected = TRUE
  )
  expect_equal(nested_overlap_description_length(nst_ndc), 56.44464320026093, tolerance = 1e-10)
  expect_equal(nested_overlap_description_length(nst_dc), 55.29410289278401, tolerance = 1e-10)
})

test_that("OverlapSBM incremental virtual_move and move deltas match exact entropy across all modes", {
  for (directed in c(FALSE, TRUE)) {
    g <- graph_create(n = 8L, directed = directed)
    edges <- matrix(c(
      1L, 2L,  1L, 2L,  1L, 3L,  2L, 3L,  2L, 2L,  3L, 4L,
      4L, 5L,  4L, 6L,  5L, 6L,  6L, 6L,  6L, 6L,  3L, 7L
    ), ncol = 2L, byrow = TRUE)
    graph_add_edges(g, edges)
    he <- matrix(c(
      1L, 1L,  1L, 2L,  1L, 1L,  1L, 2L,  1L, 1L,  2L, 3L,
      3L, 3L,  3L, 2L,  2L, 2L,  2L, 2L,  2L, 3L,  1L, 3L
    ), ncol = 2L, byrow = TRUE)

    for (dc in c(FALSE, TRUE)) {
      for (dist in c(TRUE, FALSE)) {
        if (!dc && !dist) next
        for (gen_ov in c(TRUE, FALSE)) {
          set.seed(123)
          errs <- cpp_overlap_check_moves(g$pointer, he, dc, dist, gen_ov, 80L, 1e6)
          expect_lt(max(abs(errs)), 1e-9)
        }
      }
    }
  }
})

test_that("overlap_groups and overlap_majority_groups return consistent mixed memberships", {
  g <- graph_create(n = 6L, directed = FALSE)
  edges <- matrix(c(1L, 2L,  1L, 3L,  2L, 3L,  3L, 4L,  3L, 5L,  4L, 5L),
                  ncol = 2L, byrow = TRUE)
  graph_add_edges(g, edges)
  # Vertex 3 participates in both group 1 and group 2; vertex 6 is isolated
  he <- matrix(c(1L, 1L,  1L, 1L,  1L, 1L,  2L, 2L,  2L, 2L,  2L, 2L),
               ncol = 2L, byrow = TRUE)
  st <- overlap_block_state(g, half_edges = he)
  ov <- overlap_groups(st)

  expect_equal(ov$B, 2L)
  expect_equal(nrow(ov$matrix), 6L)
  expect_equal(ncol(ov$matrix), 2L)
  # Vertex 3 has 2 half-edges in group 1 and 2 in group 2
  expect_equal(ov$matrix[3L, ], c(2L, 2L), ignore_attr = TRUE)
  # Isolated vertex 6 has 0 memberships
  expect_equal(ov$matrix[6L, ], c(0L, 0L), ignore_attr = TRUE)
  expect_equal(rowSums(ov$matrix), as.numeric(graph_degree(g)), tolerance = 1e-12)
  maj <- overlap_majority_groups(st)
  expect_length(maj, 6L)
  expect_equal(maj[1:2], c(1L, 1L))
  expect_equal(maj[4:5], c(2L, 2L))
  expect_true(is.na(maj[6L]))
})

test_that("minimize_overlap_blockmodel_dl and nested overlap recover planted overlapping communities", {
  # 4 communities (8 pure core nodes each + 6 pairwise bridges of 3 nodes = 50 nodes)
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
  ed <- list()
  for (cl in groups) {
    for (a in seq_along(cl)) {
      for (b in seq_along(cl)) {
        if (a < b && stats::runif(1) < 0.65) {
          ed[[length(ed) + 1L]] <- sort(c(cl[a], cl[b]))
        }
      }
    }
  }
  g <- graph_create(n = 50L, directed = FALSE)
  graph_add_edges(g, do.call(rbind, ed))

  set.seed(42)
  st_ov <- minimize_overlap_blockmodel_dl(g, degree_corrected = FALSE)
  st_one <- overlap_block_state(g, degree_corrected = FALSE)
  expect_s3_class(st_ov, "rgraphtool_overlap_block_state")
  expect_lt(overlap_description_length(st_ov), overlap_description_length(st_one))

  # Overlap model should strongly beat non-overlapping flat NDC SBM and detect overlapping vertices
  st_flat <- minimize_blockmodel_dl(g, degree_corrected = FALSE)
  expect_lt(overlap_description_length(st_ov), block_description_length(st_flat) - 20)
  ov_mat <- overlap_groups(st_ov)$matrix
  expect_gte(sum(rowSums(ov_mat > 0L) > 1L), 10L)

  # Nested overlap minimization
  set.seed(42)
  nst_ov <- minimize_nested_overlap_blockmodel_dl(g, degree_corrected = TRUE)
  expect_s3_class(nst_ov, "rgraphtool_nested_overlap_block_state")
  expect_true(is.finite(nested_overlap_description_length(nst_ov)))
  expect_gte(length(nst_ov$upper_levels), 1L)

  # B_min / B_max bounds respected
  st_bounded <- minimize_overlap_blockmodel_dl(g, B_min = 3L, B_max = 3L)
  expect_equal(st_bounded$overlap$B, 3L)

  # Print methods
  expect_output(print(st_ov), "overlap block state")
  expect_output(print(nst_ov), "nested overlap block state")
})

test_that("overlap_mcmc_sweep collects valid samples and entropies", {
  g <- graph_create(n = 10L, directed = FALSE)
  edges <- list()
  for (i in 1:4) for (j in (i + 1):5) edges[[length(edges) + 1L]] <- c(i, j)
  for (i in 6:9) for (j in (i + 1):10) edges[[length(edges) + 1L]] <- c(i, j)
  edges[[length(edges) + 1L]] <- c(5L, 6L)
  graph_add_edges(g, do.call(rbind, edges))

  set.seed(42)
  st <- minimize_overlap_blockmodel_dl(g, degree_corrected = TRUE)
  res <- overlap_mcmc_sweep(st, beta = 1, niter = 20L, thin = 2L, collect = TRUE)
  expect_s3_class(res$state, "rgraphtool_overlap_block_state")
  expect_length(res$majority_samples, 10L)
  expect_length(res$entropies, 10L)
  expect_true(all(is.finite(res$entropies)))
  expect_equal(res$entropies[10L], overlap_description_length(res$state), tolerance = 1e-9)

  # Greedy sweep (beta = Inf) never increases entropy
  dl0 <- overlap_description_length(st)
  res_greedy <- overlap_mcmc_sweep(st, beta = Inf, niter = 10L)
  expect_lte(overlap_description_length(res_greedy$state), dl0 + 1e-9)
})

test_that("overlap functions handle edge cases and validate inputs", {
  # Empty graph (0 edges)
  g0 <- graph_create(n = 3L, directed = FALSE)
  st0 <- overlap_block_state(g0)
  expect_equal(overlap_description_length(st0), 0)
  st0_min <- minimize_overlap_blockmodel_dl(g0)
  expect_equal(overlap_description_length(st0_min), 0)
  nst0_min <- minimize_nested_overlap_blockmodel_dl(g0)
  expect_equal(nested_overlap_description_length(nst0_min), 0)

  # Single self-loop
  g1 <- graph_create(n = 1L, directed = FALSE)
  graph_add_edges(g1, matrix(c(1L, 1L), ncol = 2L))
  st1 <- overlap_block_state(g1)
  expect_true(is.finite(overlap_description_length(st1)))

  # Invalid inputs
  g2 <- graph_create(n = 3L, directed = FALSE)
  graph_add_edges(g2, matrix(c(1L, 2L, 2L, 3L), ncol = 2L, byrow = TRUE))
  expect_error(overlap_block_state(g2, half_edges = matrix(1L, nrow = 1L, ncol = 2L)), "half_edges")
  expect_error(overlap_block_state(g2, half_edges = matrix(c(0L, 1L, 1L, 1L), ncol = 2L)), "half_edges")
  expect_error(overlap_block_state(g2, degree_corrected = NA), "degree_corrected")
  expect_error(overlap_block_state(g2, gen_overlap = NA), "gen_overlap")
  expect_error(minimize_overlap_blockmodel_dl(g2, B_min = 3L, B_max = 2L), "B_min")

  # Graph mutation invalidation
  st2 <- overlap_block_state(g2)
  graph_add_edges(g2, matrix(c(1L, 3L), ncol = 2L))
  expect_error(overlap_description_length(st2), "changed")
})
