test_that("flat model comparison returns both fits and the MDL winner", {
  g <- graph_create(12, directed = FALSE)
  graph_add_edges(g, do.call(rbind, lapply(
    list(1:6, 7:12), function(v) t(combn(v, 2))
  )))
  set.seed(43)
  result <- compare_blockmodels(g, restarts = 2, max_sweeps = 4)
  expect_s3_class(result, "graphtool_blockmodel_comparison")
  expect_equal(result$comparison$model,
               c("non-degree-corrected", "degree-corrected"))
  expect_equal(result$best_model,
               result$comparison$model[which.min(result$comparison$description_length)])
  expect_equal(
    result$comparison$description_length,
    unname(vapply(result$states, block_description_length, numeric(1))),
    tolerance = 1e-10
  )
  expect_true(result$best_state$degree_corrected ==
                identical(result$best_model, "degree-corrected"))
})

test_that("flat model comparison validates search settings", {
  expect_error(compare_blockmodels(graph_create(2), restarts = 0), "positive")
  expect_error(compare_blockmodels(graph_create(2), degree_prior = "other"), "arg")
})
