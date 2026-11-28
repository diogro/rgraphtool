
test_that("parallel_starts gives identical results for any cores and reattaches graphs", {
  skip_on_os("windows")
  g <- graph_create(40, directed = FALSE)
  set.seed(1)
  graph_add_edges(g, cbind(sample(40, 120, TRUE), sample(40, 120, TRUE)))
  set.seed(9)
  a <- parallel_starts(4, function(run) minimize_blockmodel_dl(g), cores = 1L, graph = g)
  after_a <- runif(1)
  set.seed(9)
  b <- parallel_starts(4, function(run) minimize_blockmodel_dl(g), cores = 2L, graph = g)
  after_b <- runif(1)
  expect_identical(lapply(a, `[[`, "blocks"), lapply(b, `[[`, "blocks"))
  expect_identical(after_a, after_b)  # caller's RNG stream advances the same way
  expect_equal(vapply(b, block_description_length, numeric(1)),
               vapply(a, block_description_length, numeric(1)))
  expect_error(parallel_starts(2, function(run) stop("boom"), cores = 2L), "boom")
  expect_error(parallel_starts(2, identity, cores = 0), "positive integer")
})
