deterministic_graph <- function(directed) {
  i <- 0:5999
  g <- graph_create(500, directed)
  graph_add_edges(g, cbind((i * 37) %% 500, (i * 91 + i %/% 7) %% 500) + 1)
  g
}

test_that("native description lengths match graph-tool 3.9 on a larger multigraph", {
  # Values from tools/reference/entropy_cases.py (graph-tool 3.9). The graph has
  # loops and parallel edges; the single-group undirected case exceeds the exact
  # integer-partition table and exercises the asymptotic approximation.
  expected <- list(
    "FALSE" = list(b7 = c(27048.496594042845, 27072.852331269158, 27848.159893970857),
                   b1 = c(25991.747236254370, 25661.044645038299, 26801.056938900543)),
    "TRUE" = list(b7 = c(31209.882382833843, 30585.907881270374, 32484.836649346966),
                  b1 = c(30112.022541507242, 28948.702463654743, 31402.346837379966))
  )
  old <- options(graphtool.q_cache_max = 10000)
  on.exit(options(old))
  for (directed in c(FALSE, TRUE)) {
    g <- deterministic_graph(directed)
    partitions <- list(b7 = (0:499) %% 7 + 1, b1 = rep(1, 500))
    for (name in names(partitions)) {
      b <- partitions[[name]]
      actual <- c(block_description_length(block_state(g, b, FALSE)),
                  block_description_length(block_state(g, b, TRUE)),
                  block_description_length(block_state(g, b, TRUE), "uniform"))
      expect_equal(actual, expected[[as.character(directed)]][[name]], tolerance = 1e-10)
    }
  }
})

test_that("native and R description lengths agree on random multigraphs", {
  set.seed(3)
  for (directed in c(FALSE, TRUE)) {
    g <- graph_create(30, directed)
    graph_add_edges(g, cbind(sample(30, 120, TRUE), sample(30, 120, TRUE)))
    b <- sample(c(2, 5, 9, 11), 30, TRUE)
    for (dc in c(FALSE, TRUE)) for (prior in c("distributed", "uniform")) {
      state <- block_state(g, b, dc)
      expect_equal(block_description_length(state, prior),
                   graphtool:::block_description_length_impl(state, prior),
                   tolerance = 1e-10)
    }
  }
})

test_that("integer-partition approximation matches graph-tool's log_q_approx", {
  # graph_tool.inference.libgraph_tool_inference.log_q_approx(n, k), graph-tool 3.9
  cases <- rbind(c(12000, 500, 269.42500134819056), c(30000, 3000, 432.0437406347449),
                 c(12000, 5, 29.604268704426026), c(20000, 20, 106.59874044010199),
                 c(50000, 100, 353.11435660817136), c(13188, 4941, 283.1508623375297),
                 c(13188, 13188, 283.1508623375297), c(50000, 40000, 560.8183421646834))
  for (i in seq_len(nrow(cases))) {
    expect_equal(graphtool:::cpp_log_q(cases[i, 1], cases[i, 2], 100), cases[i, 3],
                 tolerance = 1e-10)
  }
})

test_that("incremental move costs equal recomputed description lengths", {
  set.seed(4)
  for (directed in c(FALSE, TRUE)) {
    g <- graph_create(25, directed)
    graph_add_edges(g, rbind(cbind(sample(25, 90, TRUE), sample(25, 90, TRUE)),
                             c(1, 1), c(1, 1), c(2, 3), c(2, 3)))
    for (dc in c(FALSE, TRUE)) for (dist in c(FALSE, TRUE)) {
      err <- graphtool:::cpp_sbm_check_moves(graphtool:::graph_pointer(g),
                                             sample(4, 25, TRUE), dc, dist, 400L, 5000)
      expect_lt(max(abs(err)), 1e-8)
    }
  }
})

test_that("multilevel fits are reproducible and recover planted groups", {
  set.seed(9)
  g <- graph_create(60, directed = FALSE)
  planted <- rep(1:3, each = 20)
  edges <- do.call(rbind, lapply(1:400, function(i) {
    u <- sample(60, 1)
    w <- if (runif(1) < 0.9) sample(which(planted == planted[u]), 1) else sample(60, 1)
    c(u, w)
  }))
  graph_add_edges(g, edges)
  set.seed(1)
  fit <- minimize_blockmodel_dl(g)
  set.seed(1)
  again <- minimize_blockmodel_dl(g)
  expect_identical(block_membership(fit), block_membership(again))
  expect_equal(fit$fit$description_length, block_description_length(fit), tolerance = 1e-8)
  expect_lte(fit$fit$description_length, fit$fit$multilevel_description_length + 1e-8)
  expect_lte(block_description_length(fit),
             block_description_length(block_state(g, planted)) + 1e-8)
  expect_equal(length(block_sizes(fit)), 3)
  expect_equal(max(table(paste(planted, block_membership(fit)))), 20)
  expect_true(all(c("B", "entropy") %in% names(fit$fit$trace)))
})

test_that("group-count bounds are respected", {
  g <- graph_create(12, directed = FALSE)
  graph_add_edges(g, do.call(rbind, lapply(list(1:4, 5:8, 9:12), function(v) t(combn(v, 2)))))
  set.seed(2)
  expect_equal(length(block_sizes(minimize_blockmodel_dl(g, B_min = 4))), 4)
  set.seed(2)
  expect_equal(length(block_sizes(minimize_blockmodel_dl(g, B_min = 2, B_max = 2, refine = FALSE))), 2)
  expect_error(minimize_blockmodel_dl(g, B_min = 3, B_max = 2), "B_min")
  expect_error(minimize_blockmodel_dl(g, control = list(r = 1)), "between 0 and 1")
  expect_error(minimize_blockmodel_dl(g, control = list(sweeps = 1)), "named list")
})

test_that("tiny and edgeless graphs fit without error", {
  for (directed in c(FALSE, TRUE)) {
    for (n in 0:2) {
      fit <- minimize_blockmodel_dl(graph_create(n, directed))
      expect_length(block_membership(fit), n)
    }
    g <- graph_create(5, directed)
    fit <- minimize_blockmodel_dl(g)
    expect_equal(block_description_length(fit),
                 block_description_length(block_state(g)), tolerance = 1e-10)
  }
})

test_that("states from grown graphs fall back to their stored snapshot", {
  g <- graph_create(3, directed = FALSE)
  graph_add_edges(g, rbind(c(1, 2), c(2, 3)))
  state <- block_state(g, c(1, 1, 2))
  before <- block_description_length(state)
  graph_add_edges(g, rbind(c(1, 3)))
  expect_equal(block_description_length(state), before, tolerance = 1e-10)
})

test_that("merge_sweep does not stall on disconnected high-degree components", {
  g <- graph_create(120, directed = FALSE)
  edges <- do.call(rbind, lapply(0:3, function(k) {
    v <- k * 30L + seq_len(30L)
    t(combn(v, 2L))
  }))
  graph_add_edges(g, edges)
  set.seed(1)
  fit <- minimize_blockmodel_dl(g, B_min = 2, B_max = 2, control = list(c = 0.01))
  expect_equal(length(block_sizes(fit)), 2L)
  expect_true(is.finite(block_description_length(fit)))
})
