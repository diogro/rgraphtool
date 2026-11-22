test_that("flat model comparison returns both fits and the MDL winner", {
  g <- graph_create(12, directed = FALSE)
  graph_add_edges(g, do.call(rbind, lapply(
    list(1:6, 7:12), function(v) t(combn(v, 2))
  )))
  set.seed(43)
  result <- compare_blockmodels(g)
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
  expect_error(compare_blockmodels(graph_create(2), B_min = 0), "B_min")
  expect_error(compare_blockmodels(graph_create(2), degree_prior = "other"), "arg")
})

test_that("nested comparisons report log posterior odds", {
  g <- graph_create(12, directed = FALSE)
  graph_add_edges(g, do.call(rbind, lapply(list(1:4, 5:8, 9:12), function(v) t(combn(v, 2)))))
  set.seed(4)
  result <- compare_blockmodels(g, nested = TRUE, nruns = 2)
  expect_true(result$nested)
  expect_s3_class(result$states[[1]], "graphtool_nested_block_state")
  dl <- vapply(result$states, nested_description_length, numeric(1))
  expect_equal(result$comparison$description_length, unname(dl))
  expect_equal(result$log_posterior_odds, unname(dl[2] - dl[1]))
  expect_output(print(result), "ln Lambda")
  expect_error(compare_blockmodels(g, nruns = 0), "positive")
  expect_error(compare_blockmodels(g, nested = NA), "nested")
})
