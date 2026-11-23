test_that("thread configuration API validates inputs and respects options", {
  old_t <- rgraphtool_num_threads()
  old_opt <- getOption("rgraphtool.nthreads", NULL)
  on.exit({
    options(rgraphtool.nthreads = old_opt)
    rgraphtool_set_num_threads(old_t)
  }, add = TRUE)

  hw <- rgraphtool_max_threads()
  expect_true(is.integer(hw) && length(hw) == 1L && hw >= 1L)

  rgraphtool_set_num_threads(2L)
  expect_equal(rgraphtool_num_threads(), 2L)
  rgraphtool_set_num_threads("auto")
  expect_equal(rgraphtool_num_threads(), hw)
  rgraphtool_set_num_threads(1L)
  expect_equal(rgraphtool_num_threads(), 1L)

  options(rgraphtool.nthreads = 3L)
  expect_equal(rgraphtool_num_threads(), 3L)
  options(rgraphtool.nthreads = "auto")
  expect_equal(rgraphtool_num_threads(), hw)
  options(rgraphtool.nthreads = NULL)
  expect_equal(rgraphtool_num_threads(), 1L)

  expect_error(rgraphtool_set_num_threads(0L), "positive integer")
  expect_error(rgraphtool_set_num_threads(-2L), "positive integer")
  expect_error(rgraphtool_set_num_threads(1.5), "positive integer")
  expect_error(rgraphtool_set_num_threads("invalid"), "positive integer")
})

test_that("strictly agglomerative search, simple_init, warm-start, and parallel threads work", {
  make_modular <- function(n, B, k_avg = 10, directed = FALSE) {
    sz <- n %/% B
    m_in <- as.integer(round(n * k_avg * 0.42))
    m_out <- as.integer(round(n * k_avg * 0.08))
    g_id <- sample.int(B, m_in, replace = TRUE)
    u_in <- (g_id - 1L) * sz + sample.int(sz, m_in, replace = TRUE)
    v_in <- (g_id - 1L) * sz + sample.int(sz, m_in, replace = TRUE)
    keep <- u_in != v_in & u_in <= n & v_in <= n
    u_out <- sample.int(n, m_out, replace = TRUE)
    v_out <- sample.int(n, m_out, replace = TRUE)
    keep2 <- u_out != v_out
    ed <- rbind(cbind(u_in[keep], v_in[keep]), cbind(u_out[keep2], v_out[keep2]))
    g <- graph_add_edges(graph_create(n, directed = directed), ed)
    truth <- rep(seq_len(B), each = sz)[seq_len(n)]
    list(g = g, truth = truth)
  }

  set.seed(11L)
  dat <- make_modular(2400L, 12L, k_avg = 10, directed = FALSE)
  g <- dat$g
  truth <- dat$truth

  # 1. Strictly agglomerative flat search (bisection = FALSE) with 1 vs 2 threads
  set.seed(21L)
  st_a1 <- minimize_blockmodel_dl(g, nthreads = 1L,
                                  control = list(bisection = FALSE, r = 0.5))
  expect_equal(st_a1$fit$description_length, block_description_length(st_a1), tolerance = 1e-9)
  expect_gt(partition_overlap(block_membership(st_a1), truth, norm = TRUE), 0.95)

  set.seed(21L)
  st_a2 <- minimize_blockmodel_dl(g, nthreads = 2L,
                                  control = list(bisection = FALSE, r = 0.5))
  expect_equal(st_a2$fit$description_length, block_description_length(st_a2), tolerance = 1e-9)
  expect_gt(partition_overlap(block_membership(st_a2), truth, norm = TRUE), 0.95)

  # Parallel determinism: same seed + same nthreads gives identical result
  set.seed(21L)
  st_a2b <- minimize_blockmodel_dl(g, nthreads = 2L,
                                   control = list(bisection = FALSE, r = 0.5))
  expect_equal(st_a2$fit$description_length, st_a2b$fit$description_length, tolerance = 1e-12)
  expect_equal(block_membership(st_a2), block_membership(st_a2b))

  # 2. Warm-start flat minimization from integer vector and from block state
  noisy <- truth
  noisy[seq_len(120L)] <- sample.int(12L, 120L, replace = TRUE)
  st_ws <- minimize_blockmodel_dl(g, init_blocks = noisy, nthreads = 2L)
  expect_lte(block_description_length(st_ws),
             block_description_length(block_state(g, noisy)))
  expect_gt(partition_overlap(block_membership(st_ws), truth, norm = TRUE), 0.95)
  st_ws2 <- minimize_blockmodel_dl(g, init_blocks = st_ws, refine = FALSE)
  expect_equal(block_description_length(st_ws2), block_description_length(st_ws),
               tolerance = 1e-9)

  # 3. Nested simple_init (bottom-up and top-down) and warm-start init_levels
  set.seed(31L)
  nst_si <- minimize_nested_blockmodel_dl(
    g, simple_init = TRUE, base_factor = 10, B_min_base = 60L,
    top_factor = 2, strict = FALSE, nthreads = 2L, control = list(r = 0.5)
  )
  expect_equal(nst_si$fit$description_length, nested_description_length(nst_si),
               tolerance = 1e-8)
  expect_gt(partition_overlap(nested_levels(nst_si)[[1L]], truth, norm = TRUE), 0.95)

  set.seed(32L)
  nst_rev <- minimize_nested_blockmodel_dl(
    g, simple_init = TRUE, simple_init_reverse = TRUE,
    base_factor = 10, B_min_base = 40L, top_factor = 2,
    strict = TRUE, refine = TRUE, nthreads = 2L, control = list(r = 0.5)
  )
  expect_equal(nst_rev$fit$description_length, nested_description_length(nst_rev),
               tolerance = 1e-8)
  expect_gt(partition_overlap(nested_levels(nst_rev)[[1L]], truth, norm = TRUE), 0.90)

  nst_ws <- minimize_nested_blockmodel_dl(g, init_levels = nst_si, nthreads = 2L)
  expect_lte(nested_description_length(nst_ws), nested_description_length(nst_si) + 1e-7)

  # 4. Directed graph with parallel threads + simple_init + bisection = FALSE
  set.seed(41L)
  dat_d <- make_modular(400L, 4L, k_avg = 12, directed = TRUE)
  st_d <- minimize_blockmodel_dl(dat_d$g, nthreads = 2L,
                                 control = list(bisection = FALSE, r = 0.5))
  expect_equal(st_d$fit$description_length, block_description_length(st_d), tolerance = 1e-9)
  expect_gt(partition_overlap(block_membership(st_d), dat_d$truth, norm = TRUE), 0.90)
  nst_d <- minimize_nested_blockmodel_dl(dat_d$g, simple_init = TRUE, B_min_base = 20L,
                                         nthreads = 2L, control = list(r = 0.5))
  expect_equal(nst_d$fit$description_length, nested_description_length(nst_d), tolerance = 1e-8)
})

test_that("NestedOverlap incremental upper-level updates match full recomputation", {
  set.seed(99L)
  edges <- rbind(
    c(1L, 2L), c(1L, 2L), c(2L, 3L), c(3L, 4L), c(4L, 5L),
    c(5L, 6L), c(6L, 1L), c(1L, 1L), c(3L, 6L), c(2L, 5L),
    c(4L, 7L), c(7L, 8L), c(7L, 8L), c(8L, 1L)
  )
  for (dir in c(FALSE, TRUE)) {
    g <- graph_add_edges(graph_create(8L, directed = dir), edges)
    nst <- minimize_nested_overlap_blockmodel_dl(g, degree_corrected = TRUE)
    errs <- rgraphtool:::cpp_nested_overlap_check_moves(
      rgraphtool:::graph_pointer(g), nst$half_edges, nst$upper_levels,
      nst$degree_corrected, TRUE, TRUE, 300L, 5000
    )
    expect_lt(max(errs), 1e-10)
  }
})

test_that("large-network controls validate edge cases and invalid inputs", {
  g0 <- graph_create(0L, directed = FALSE)
  expect_equal(length(block_membership(
    minimize_blockmodel_dl(g0, nthreads = 2L, control = list(bisection = FALSE))
  )), 0L)

  g1 <- graph_create(1L, directed = FALSE)
  expect_equal(block_membership(
    minimize_blockmodel_dl(g1, nthreads = 2L, control = list(bisection = FALSE))
  ), 1L)
  expect_equal(nested_levels(
    minimize_nested_blockmodel_dl(g1, simple_init = TRUE, nthreads = 2L)
  )[[1L]], 1L)

  g_small <- graph_add_edges(graph_create(6L, directed = FALSE),
                             rbind(c(1L, 1L), c(1L, 2L), c(1L, 2L), c(3L, 4L), c(5L, 6L)))
  expect_error(minimize_blockmodel_dl(g_small, control = list(bisection = NA)),
               "TRUE or FALSE")
  expect_error(minimize_blockmodel_dl(g_small, nthreads = 0L), "positive integer")
  expect_error(minimize_blockmodel_dl(g_small, init_blocks = c(1L, 2L)), "length")
  expect_error(minimize_nested_blockmodel_dl(g_small, simple_init = NA), "TRUE or FALSE")
  expect_error(minimize_nested_blockmodel_dl(g_small, base_factor = 0), "positive finite")
  expect_error(minimize_nested_blockmodel_dl(g_small, B_min_base = 0L), "positive integer")
  expect_error(minimize_nested_blockmodel_dl(g_small, top_factor = 1), "greater than 1")
  expect_error(minimize_nested_blockmodel_dl(g_small, strict = NA), "TRUE or FALSE")
})
