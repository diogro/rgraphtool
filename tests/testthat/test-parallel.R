
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

test_that("forked children can run threaded code after the parent's pool started", {
  skip_on_os("windows")
  skip_on_cran()
  set.seed(1)
  f <- matrix(rnorm(400 * 5), 400, 5)
  x <- f %*% matrix(rnorm(5 * 1500), 5, 1500) + matrix(rnorm(400 * 1500), 400, 1500)
  # parent spawns pool workers; a child inherits the pool object but not its threads
  p <- shrink_pcor(x, top_k = 500, nthreads = 3)
  res <- parallel::mclapply(1:2, function(i) {
    identical(shrink_pcor(x, top_k = 500, nthreads = 3)$top, p$top)
  }, mc.cores = 2)
  expect_identical(unlist(res), c(TRUE, TRUE))
  expect_identical(shrink_pcor(x, top_k = 500, nthreads = 3)$top, p$top)
})
