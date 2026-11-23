test_that("description lengths match graph-tool 3.9 reference cases", {
  edges <- rbind(c(1, 2), c(1, 2), c(2, 1), c(3, 3), c(3, 2))
  expected <- list(
    c(12.86636600121904, 12.108680299521524, 12.48597453066299),
    c(14.134877326682545, 13.559513181778986, 14.475803913653138)
  )
  for (i in seq_along(expected)) {
    g <- graph_create(4, directed = i == 2)
    graph_add_edges(g, edges)
    noncorrected <- block_state(g, c(2, 1, 2, 1), FALSE)
    corrected <- block_state(g, c(2, 1, 2, 1), TRUE)
    actual <- c(block_description_length(noncorrected),
                block_description_length(corrected),
                block_description_length(corrected, "uniform"))
    expect_equal(actual, expected[[i]], tolerance = 1e-10)
  }
})

test_that("description length covers empty graphs and isolated vertices", {
  for (directed in c(FALSE, TRUE)) {
    expect_equal(block_description_length(block_state(graph_create(0, directed))), 0)
    state <- block_state(graph_create(2, directed))
    expect_equal(block_description_length(state), log(2))
    expect_equal(block_description_length(state, "uniform"), log(2))
  }
  expect_error(rgraphtool:::cpp_log_partition_table(10000L, 1000L), "10 million")
})

test_that("multilevel fitting is reproducible and improves on one group", {
  g <- graph_create(12, directed = FALSE)
  graph_add_edges(g, do.call(rbind, lapply(list(1:6, 7:12), function(v) t(combn(v, 2)))))
  before <- block_description_length(block_state(g))
  set.seed(12)
  fit <- minimize_blockmodel_dl(g)
  set.seed(12)
  repeat_fit <- minimize_blockmodel_dl(g)
  expect_equal(block_membership(fit), block_membership(repeat_fit))
  expect_lt(block_description_length(fit), before)
  expect_equal(sort(unname(block_sizes(fit))), c(6, 6))
  expect_equal(graph_num_edges(g), 30)
})

test_that("search remembers its degree prior and rejects invalid settings", {
  g <- graph_create(3)
  graph_add_edges(g, rbind(c(1, 2), c(2, 3)))
  set.seed(5)
  fit <- minimize_blockmodel_dl(g, degree_prior = "uniform")
  expect_equal(block_description_length(fit), block_description_length(fit, "uniform"))
  expect_error(minimize_blockmodel_dl(g, refine = NA), "refine")
  expect_error(minimize_blockmodel_dl(g, control = list(mh_sweeps = -1)), "non-negative")
  expect_error(block_description_length(fit, "other"), "arg")
})
