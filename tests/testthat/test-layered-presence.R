# Given per-layer vertex sets (`present`) and weight-unit invariance of the
# independent-vs-edge-covariates comparison.

presence_graph <- function(directed = FALSE) {
  set.seed(71L)
  n <- 24L
  e1 <- t(utils::combn(1:20, 2))[stats::runif(190) < 0.3, ]
  e2 <- t(utils::combn(5:24, 2))[stats::runif(190) < 0.3, ]
  ed <- rbind(e1, e2)
  lay <- rep(c("a", "b"), c(nrow(e1), nrow(e2)))
  w <- stats::rnorm(nrow(ed), 0.5, 1)
  g <- graph_add_edges(graph_create(n, directed), ed,
                       properties = list(layer = lay, w = w, we = exp(w)))
  pres <- cbind(a = seq_len(n) <= 22L, b = seq_len(n) >= 3L)
  list(g = g, present = pres, blocks = rep(1:4, each = 6L))
}

test_that("given presence sums to one over all graphs when every vertex is present", {
  N <- 3L
  pairs <- list()
  for (i in 1:N) for (j in i:N) pairs[[length(pairs) + 1L]] <- c(i, j)
  b <- c(1L, 1L, 2L)
  n <- tabulate(b)
  # the partition prior is not part of the graph likelihood
  pdl <- log(N) + lchoose(N - 1, 1) + lfactorial(N) - sum(lfactorial(n))
  for (dc in c(TRUE, FALSE)) {
    tot <- c(layered = 0, pp = 0)
    for (p1 in seq_along(pairs)) for (p2 in seq_along(pairs)) {
      g <- graph_add_edges(graph_create(N, FALSE), rbind(pairs[[p1]], pairs[[p2]]),
                           properties = list(layer = 1:2))
      pres <- matrix(TRUE, N, 2L)
      s_l <- layered_description_length(layered_block_state(
        g, "layer", TRUE, blocks = b, degree_corrected = dc, present = pres))
      s_p <- pp_description_length(pp_block_state(
        g, b, degree_corrected = dc, layers = "layer", independent = TRUE, present = pres))
      tot <- tot + exp(-(c(s_l, s_p) - pdl))
    }
    expect_equal(unname(tot), c(1, 1), tolerance = 1e-10)
  }
})

test_that("presence is validated against the edges and the layer levels", {
  x <- presence_graph()
  bad <- x$present
  bad[10L, "b"] <- FALSE  # vertex 10 has edges in layer b
  expect_error(layered_block_state(x$g, "layer", TRUE, blocks = x$blocks, present = bad),
               "not marked present")
  expect_error(pp_block_state(x$g, x$blocks, layers = "layer", independent = TRUE,
                              present = bad), "not marked present")
  expect_error(layered_block_state(x$g, "layer", TRUE, blocks = x$blocks,
                                   present = x$present[, 1, drop = FALSE]), "logical matrix")
  renamed <- x$present
  colnames(renamed) <- c("a", "c")
  expect_error(layered_block_state(x$g, "layer", TRUE, blocks = x$blocks, present = renamed),
               "layer levels")
  # named columns are matched to layers, not taken in order
  swapped <- x$present[, c("b", "a")]
  expect_equal(
    layered_description_length(layered_block_state(x$g, "layer", TRUE, blocks = x$blocks,
                                                   present = swapped)),
    layered_description_length(layered_block_state(x$g, "layer", TRUE, blocks = x$blocks,
                                                   present = x$present)))
})

test_that("presence-conditioned moves match full recomputation", {
  for (directed in c(FALSE, TRUE)) {
    x <- presence_graph(directed)
    st <- layered_block_state(x$g, "layer", TRUE, blocks = x$blocks, present = x$present)
    for (dc in c(TRUE, FALSE)) {
      set.seed(3L)
      err <- rgraphtool:::cpp_layered_check_moves(
        rgraphtool:::graph_pointer(x$g), st$blocks, st$layers, TRUE, dc, TRUE, 40L, 1000,
        present = x$present)
      expect_lt(max(abs(err)), 1e-8)
      set.seed(3L)
      err <- rgraphtool:::cpp_pp_check_moves(
        rgraphtool:::graph_pointer(x$g), x$blocks, FALSE, dc, TRUE, 40L, 1000,
        layers = st$layers, independent = TRUE, present = x$present)
      expect_lt(max(abs(err)), 1e-8)
    }
  }
})

test_that("presence-conditioned fits run and keep present", {
  x <- presence_graph()
  set.seed(5L)
  fit <- minimize_layered_blockmodel_dl(x$g, "layer", TRUE, present = x$present)
  expect_identical(fit$present, unname(x$present))
  expect_equal(layered_description_length(fit),
               layered_description_length(layered_block_state(
                 x$g, "layer", TRUE, blocks = block_membership(fit), present = x$present)))
})

test_that("independent minus edge-covariates DL does not depend on weight units", {
  x <- presence_graph()
  diff_dl <- function(scale, w, tp, model) {
    g <- x$g
    graph_set_edge_property(g, "ws", scale * graph_edge_property(g, w))
    f <- function(ind) switch(model,
      layered = layered_description_length(layered_block_state(
        g, "layer", ind, blocks = x$blocks, weights = "ws", types = tp)),
      pp = pp_description_length(pp_block_state(
        g, x$blocks, layers = "layer", independent = ind, weights = "ws", types = tp)))
    f(TRUE) - f(FALSE)
  }
  for (tp in c("real-normal", "real-exponential")) {
    w <- if (tp == "real-normal") "w" else "we"
    for (model in c("layered", "pp")) {
      d <- vapply(c(1, 10, 1000), diff_dl, numeric(1), w = w, tp = tp, model = model)
      expect_equal(d, rep(d[1], 3), tolerance = 1e-8)
    }
  }
})

test_that("compare_layered_models needs two layers", {
  g <- graph_add_edges(graph_create(4, FALSE), rbind(c(1, 2), c(3, 4)),
                       properties = list(layer = c(1L, 1L)))
  expect_error(compare_layered_models(g, "layer"), "at least 2 layers")
})

test_that("compare_layered_models accepts presence with layer-name columns", {
  x <- presence_graph()
  set.seed(9L)
  cmp <- compare_layered_models(x$g, "layer", nested = FALSE, present = x$present)
  ind <- cmp$states[["independent_layers"]]
  expect_identical(ind$present, unname(x$present))
  expect_true(all(is.finite(cmp$summary$description_length)))
})
