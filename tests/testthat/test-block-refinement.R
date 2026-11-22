test_that("refinement never increases description length", {
  g <- graph_create(12, directed = FALSE)
  graph_add_edges(g, do.call(rbind, lapply(
    list(1:6, 7:12), function(v) t(combn(v, 2))
  )))
  initial <- block_state(g, rep(1:2, each = 6))
  set.seed(18)
  refined <- refine_blockmodel_dl(g, initial)
  expect_lte(block_description_length(refined),
             block_description_length(initial) + 1e-10)
  expect_equal(graph_num_edges(g), 30)
})

test_that("refinement can improve a poor initial partition", {
  g <- graph_create(12, directed = FALSE)
  graph_add_edges(g, do.call(rbind, lapply(
    list(1:6, 7:12), function(v) t(combn(v, 2))
  )))
  initial <- block_state(g, c(1, 2, 1, 1, 2, 1, 1, 1, 2, 2, 1, 1))
  set.seed(31)
  refined <- refine_blockmodel_dl(g, initial)
  expect_lt(block_description_length(refined), block_description_length(initial))
})

test_that("refinement validates the graph, state, and settings", {
  g <- graph_create(2)
  expect_error(refine_blockmodel_dl(g, block_state(graph_create(3))),
               "same number of vertices")
  expect_error(refine_blockmodel_dl(g, list(blocks = 1:2)), "block_state")
  expect_error(refine_blockmodel_dl(g, block_state(g), epsilon = -1), "non-negative")
})
