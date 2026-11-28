deterministic_graph <- function(directed) {
  i <- 0:5999
  g <- graph_create(500, directed)
  graph_add_edges(g, cbind((i * 37) %% 500, (i * 91 + i %/% 7) %% 500) + 1)
  g
}

hierarchy7 <- list((0:499) %% 7 + 1, c(1, 1, 2, 2, 3, 3, 3), c(1, 1, 2), c(1, 1), 1)

test_that("nested description lengths match graph-tool 3.9", {
  # tools/reference/nested_entropy_cases.py with hierarchy7 (graph-tool 3.9)
  expected <- list(
    "FALSE" = c(27066.766756086003, 27091.122493312316, 27866.430056014018),
    "TRUE" = c(31231.814072219378, 30607.839570655913, 32506.768338732512)
  )
  old <- options(rgraphtool.q_cache_max = 10000)
  on.exit(options(old))
  for (directed in c(FALSE, TRUE)) {
    g <- deterministic_graph(directed)
    ndc <- nested_block_state(g, hierarchy7, degree_corrected = FALSE)
    dc <- nested_block_state(g, hierarchy7)
    actual <- c(nested_description_length(ndc), nested_description_length(dc),
                nested_description_length(dc, "uniform"))
    expect_equal(actual, expected[[as.character(directed)]], tolerance = 1e-10)
  }
})

test_that("a one-level hierarchy equals the flat description length", {
  g <- deterministic_graph(FALSE)
  for (dc in c(FALSE, TRUE)) {
    flat <- block_description_length(block_state(g, hierarchy7[[1]], dc))
    nested <- nested_description_length(nested_block_state(g, hierarchy7[1], dc))
    expect_equal(nested, flat, tolerance = 1e-10)
  }
})

test_that("incremental moves at every level equal recomputed nested scores", {
  set.seed(5)
  for (directed in c(FALSE, TRUE)) {
    g <- graph_create(60, directed)
    graph_add_edges(g, rbind(cbind(sample(60, 250, TRUE), sample(60, 250, TRUE)),
                             c(1, 1), c(2, 2), c(3, 4), c(3, 4)))
    levels <- list(sample(8, 60, TRUE), c(1, 1, 2, 2, 3, 3, 1, 2), c(1, 1, 2), c(1, 1), 1)
    levels <- rgraphtool:::validate_hierarchy(levels, 60)
    for (dc in c(FALSE, TRUE)) for (dist in c(FALSE, TRUE)) {
      err <- rgraphtool:::cpp_nested_check_moves(rgraphtool:::graph_pointer(g), levels,
                                                dc, dist, 1000L, 5000)
      expect_lt(max(abs(err)), 1e-8)
    }
  }
})

test_that("hierarchies are validated, relabelled, and projected", {
  g <- graph_create(4, directed = FALSE)
  graph_add_edges(g, rbind(c(1, 2), c(3, 4)))
  state <- nested_block_state(g, list(c(5, 5, 9, 9), c(3, 3)))
  expect_equal(nested_levels(state), list(c(1L, 1L, 2L, 2L), c(1L, 1L)))
  expect_equal(nested_project(state, 1), rep(1L, 4))
  expect_equal(block_membership(nested_level_state(state, 0)), c(1L, 1L, 2L, 2L))
  expect_equal(nested_summary(state)$B, c(2L, 1L))
  expect_equal(sum(nested_level_entropy(state)), nested_description_length(state))
  expect_error(nested_block_state(g, list(c(1, 1, 2, 2), c(1, 2, 3))), "levels\\[\\[2\\]\\]")
  expect_error(nested_block_state(g, list(c(1, 1, 2))), "per vertex")
  expect_error(nested_project(state, 2), "between 0 and 1")
  expect_output(print(state), "l: 0, N: 4, B: 2")
  graph_add_edges(g, rbind(c(1, 3)))
  expect_error(nested_description_length(state), "changed")
})

test_that("nested fits are reproducible and find a planted hierarchy", {
  set.seed(11)
  n <- 96
  planted <- rep(1:8, each = 12)
  super <- (planted + 1) %/% 2
  edges <- do.call(rbind, lapply(1:1500, function(i) {
    u <- sample(n, 1)
    p <- runif(1)
    pool <- if (p < 0.75) which(planted == planted[u]) else if (p < 0.92)
      which(super == super[u]) else seq_len(n)
    c(u, pool[sample.int(length(pool), 1)])
  }))
  g <- graph_create(n, directed = FALSE)
  graph_add_edges(g, edges)
  set.seed(3)
  fit <- minimize_nested_blockmodel_dl(g)
  set.seed(3)
  again <- minimize_nested_blockmodel_dl(g)
  expect_identical(nested_levels(fit), nested_levels(again))
  expect_equal(fit$fit$description_length, nested_description_length(fit), tolerance = 1e-8)
  expect_equal(nested_summary(fit)$B[1], 8L)
  # accepted zero-temperature level searches and refinement never worsen the fit
  expect_true(all(diff(fit$fit$trace) <= 1e-8))
  truth <- nested_block_state(g, list(planted, (seq_len(8) + 1) %/% 2, rep(1, 4)))
  expect_lte(nested_description_length(fit), nested_description_length(truth) + 1e-8)
  set.seed(3)
  flat <- minimize_blockmodel_dl(g)
  expect_equal(length(block_sizes(flat)), 8)
})

test_that("tiny graphs fit and bounds apply to the lowest level", {
  for (directed in c(FALSE, TRUE)) {
    fit <- minimize_nested_blockmodel_dl(graph_create(1, directed))
    expect_equal(nested_levels(fit), list(1L))
    g <- graph_create(12, directed)
    graph_add_edges(g, do.call(rbind, lapply(list(1:4, 5:8, 9:12), function(v) t(combn(v, 2)))))
    set.seed(1)
    fit <- minimize_nested_blockmodel_dl(g, B_min = 2, B_max = 2)
    expect_equal(max(nested_levels(fit)[[1]]), 2L)
  }
  expect_error(minimize_nested_blockmodel_dl(graph_create(0)), "at least one vertex")
})

test_that("the hierarchy plot runs and returns its layout", {
  g <- graph_create(4, directed = FALSE)
  graph_add_edges(g, rbind(c(1, 2), c(3, 4), c(2, 3)))
  state <- nested_block_state(g, list(c(1, 1, 2, 2), c(1, 1)))
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  on.exit({ grDevices::dev.off(); unlink(file) })
  layout <- plot_hierarchy(state)
  expect_equal(dim(layout$vertices), c(4L, 2L))
  expect_equal(nrow(layout$groups[[1]]), 2L)
  expect_error(plot_hierarchy(state, beta = 2), "beta")
})

test_that("C. elegans nested fits reach graph-tool 3.9 description lengths", {
  skip_on_cran()
  celegans <- tryCatch(netzschleuder_network("celegansneural"), error = function(e) NULL)
  skip_if(is.null(celegans), "C. elegans network could not be downloaded")
  expect_true(graph_num_vertices(celegans$graph) == 297)
  expect_true(graph_num_edges(celegans$graph) == 2359)
  # graph-tool 3.9, seeds 0-9 (tools/reference/celegans_nested_upstream.csv):
  # non-degree-corrected 8459.65-8615.34, degree-corrected 8539.59-8662.94.
  upstream_max <- c("FALSE" = 8615.34, "TRUE" = 8662.94)
  for (dc in c(FALSE, TRUE)) {
    scores <- vapply(1:2, function(seed) {
      set.seed(seed)
      nested_description_length(minimize_nested_blockmodel_dl(celegans$graph,
                                                              degree_corrected = dc))
    }, numeric(1))
    expect_lte(min(scores), upstream_max[[as.character(dc)]])
  }
})

test_that("nested_block_state keeps parents when relabelling user hierarchies", {
  g <- graph_create(120L, directed = FALSE)
  set.seed(4)
  graph_add_edges(g, cbind(sample.int(120L, 400L, TRUE), sample.int(120L, 400L, TRUE)))
  l1 <- rep(c(3L, 1L, 4L, 2L), each = 30L)
  l2 <- c(1L, 2L, 1L, 2L)  # groups 1,3 -> parent 1; groups 2,4 -> parent 2
  st <- nested_block_state(g, list(l1, l2, c(1L, 1L)))
  # the same hierarchy in first-appearance order: 3,1,4,2 -> 1,2,3,4
  st_ref <- nested_block_state(g, list(rep(1:4, each = 30L), c(1L, 2L, 1L, 2L)[c(3L, 1L, 4L, 2L)],
                                       c(1L, 1L)))
  expect_equal(nested_description_length(st), nested_description_length(st_ref))
  expect_equal(as.integer(nested_project(st, level = 1L)),
               as.integer(nested_project(st_ref, level = 1L)))
})
