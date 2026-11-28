.ref_shrink_pcor <- function(x, lambda = NULL, center = TRUE) {
  N <- nrow(x); M <- ncol(x)
  xc <- if (center) x - rowMeans(x) else x
  df <- if (center) M - 1 else M
  xs <- xc / sqrt(rowSums(xc^2) / df)
  R <- tcrossprod(xs) / df
  if (is.null(lambda)) {
    w <- array(0, c(N, N, M))
    for (k in seq_len(M)) w[, , k] <- tcrossprod(xs[, k])
    wbar <- apply(w, c(1, 2), mean)
    ssd <- apply(sweep(w, c(1, 2), wbar)^2, c(1, 2), sum)
    v <- if (center) M / (M - 1)^3 * ssd else ssd / (M * (M - 1))
    off <- row(R) != col(R)
    lambda <- min(1, max(0, sum(v[off]) / sum(R[off]^2)))
  }
  Rs <- (1 - lambda) * R + lambda * diag(N)
  P <- solve(Rs)
  pc <- -P / sqrt(outer(diag(P), diag(P)))
  diag(pc) <- 0
  list(matrix = pc, lambda = lambda)
}

test_that("shrink_pcor matches the explicit Schafer-Strimmer estimator", {
  set.seed(3)
  for (dims in list(c(12, 30), c(25, 10))) {
    x <- matrix(rnorm(prod(dims)), dims[1], dims[2])
    x[1, ] <- x[1, ] + x[2, ]
    for (cen in c(TRUE, FALSE)) {
      ref <- .ref_shrink_pcor(x, center = cen)
      got <- shrink_pcor(x, center = cen)
      expect_equal(got$lambda, ref$lambda, tolerance = 1e-10)
      expect_equal(got$matrix, ref$matrix, ignore_attr = TRUE, tolerance = 1e-10)
      ref_f <- .ref_shrink_pcor(x, lambda = 0.25, center = cen)
      expect_equal(shrink_pcor(x, lambda = 0.25, center = cen)$matrix, ref_f$matrix,
                   ignore_attr = TRUE, tolerance = 1e-10)
    }
  }
})

test_that("shrink_pcor pairs, top_k and edge cases", {
  set.seed(4)
  x <- matrix(rnorm(40 * 15), 40, 15)
  full <- shrink_pcor(x, lambda = 0.4)$matrix
  pr <- cbind(c(1L, 3L, 10L), c(2L, 40L, 11L))
  expect_equal(shrink_pcor(x, lambda = 0.4, pairs = pr)$pairs, full[pr], tolerance = 1e-12)
  top <- shrink_pcor(x, lambda = 0.4, top_k = 30L)$top
  ut <- which(upper.tri(full))
  expect_equal(abs(top$pcor), sort(abs(full[ut]), decreasing = TRUE)[1:30], tolerance = 1e-12)
  expect_true(all(top$from < top$to))
  expect_null(shrink_pcor(x, top_k = 5L)$matrix)
  expect_equal(shrink_pcor_lambda(x), shrink_pcor(x, top_k = 1L)$lambda)
  expect_true(all(shrink_pcor(x, lambda = 1)$matrix == 0))
  y <- matrix(rnorm(5 * 50), 5, 50)
  P <- solve(cor(t(y)))
  pc0 <- -P / sqrt(outer(diag(P), diag(P)))
  diag(pc0) <- 0
  expect_equal(shrink_pcor(y, lambda = 0)$matrix, pc0, ignore_attr = TRUE, tolerance = 1e-10)
  expect_error(shrink_pcor(x, lambda = 0), "singular")
  x[2, ] <- 1
  expect_error(shrink_pcor(x), "zero or non-finite variance")
})

test_that("shrink_pcor top_k does not depend on the thread count", {
  set.seed(5)
  x <- matrix(rnorm(700 * 30), 700, 30)
  a <- shrink_pcor(x, top_k = 500L, nthreads = 1L)$top
  b <- shrink_pcor(x, top_k = 500L, nthreads = 4L)$top
  expect_identical(a, b)
})
