# Merge-split samplers against exact enumeration of all set partitions.
# The weighted, planted-partition and layered merge-split moves used wrong
# proposal probabilities and a state-dependent number of attempts per sweep
# (TV 0.07-0.15 vs exact); with the corrected moves TV is at the sampling
# noise level of single-node moves.
.set_partitions <- function(n) {
  out <- list()
  rec <- function(v) {
    if (length(v) == n) { out[[length(out) + 1L]] <<- v; return(invisible()) }
    for (k in seq_len(max(v) + 1L)) rec(c(v, k))
  }
  rec(1L)
  out
}
.canon <- function(b) paste(match(b, unique(b)), collapse = "")
.posterior_tv <- function(mk, dl, n, nsamples, seed = 7L) {
  parts <- .set_partitions(n)
  S <- vapply(parts, function(b) dl(mk(as.integer(b))), numeric(1))
  p <- exp(-(S - min(S)))
  p <- p / sum(p)
  names(p) <- vapply(parts, .canon, "")
  set.seed(seed)
  sp <- sample_posterior(mk(rep(1L, n)), nsamples = nsamples, thin = 2L, burnin = 100L,
                         multiflip = TRUE)
  lab <- vapply(sp$samples, .canon, "")
  emp <- as.numeric(table(factor(lab, levels = names(p)))) / length(lab)
  0.5 * sum(abs(emp - p))
}

test_that("weighted merge-split samples the exact posterior", {
  N <- 6L
  e <- rbind(c(1, 2), c(2, 3), c(1, 3), c(4, 5), c(5, 6), c(4, 6), c(3, 4), c(2, 2))
  g <- graph_add_edges(graph_create(N, directed = FALSE), e)
  set.seed(2)
  graph_set_edge_property(g, "w", rexp(graph_num_edges(g)))
  mk <- function(b) weighted_block_state(g, "w", "real-exponential", blocks = b,
                                         degree_corrected = FALSE)
  expect_lt(.posterior_tv(mk, weighted_description_length, N, 15000L), 0.03)
})

test_that("planted-partition and layered merge-split sample the exact posterior", {
  N <- 5L
  ed <- rbind(c(1, 2), c(2, 3), c(1, 3), c(3, 4), c(4, 5), c(2, 5))
  g <- graph_add_edges(graph_create(N, directed = FALSE), ed,
                       properties = list(layer = c(1L, 1L, 2L, 2L, 1L, 2L)))
  dl <- function(st) rgraphtool:::state_description_length(st)
  expect_lt(.posterior_tv(function(b) pp_block_state(g, b), dl, N, 15000L), 0.03)
  expect_lt(.posterior_tv(function(b) layered_block_state(g, "layer", independent = FALSE, blocks = b),
                          dl, N, 15000L), 0.03)
  expect_lt(.posterior_tv(function(b) layered_block_state(g, "layer", independent = TRUE, blocks = b),
                          dl, N, 15000L), 0.03)
})

test_that("nested weighted MCMC changes B and samples the exact hierarchy posterior", {
  # Regression: level moves never created or removed groups, so B was fixed
  # at every level and samples were conditional on the starting hierarchy.
  N <- 5L
  e <- rbind(c(1, 2), c(2, 3), c(1, 3), c(3, 4), c(4, 5), c(2, 5), c(1, 4))
  g <- graph_add_edges(graph_create(N, directed = FALSE), e)
  set.seed(3)
  w <- rexp(nrow(e))
  mk <- function(b0, b1) nested_weighted_block_state(
    g, weights = w, types = "real-exponential",
    levels = list(b0, b1, rep(1L, max(b1))), degree_corrected = FALSE
  )
  key <- function(b0, b1) {
    u <- unique(b0)
    b1r <- b1[u]
    paste(.canon(b0), .canon(b1r), sep = "|")
  }
  H <- list()
  for (b0 in .set_partitions(N)) for (b1 in .set_partitions(max(b0))) H[[length(H) + 1L]] <- list(b0, b1)
  S <- vapply(H, function(h) nested_weighted_description_length(mk(h[[1]], h[[2]])), numeric(1))
  p <- exp(-(S - min(S)))
  p <- p / sum(p)
  names(p) <- vapply(H, function(h) key(h[[1]], h[[2]]), "")
  set.seed(11)
  sp <- sample_posterior(mk(rep(1L, N), 1L), nsamples = 12000L, thin = 1L, burnin = 200L,
                         mcmc_args = list(d = 0.3))
  ks <- vapply(sp$samples, function(s) key(s[[1]], s[[2]]), "")
  expect_true(all(ks %in% names(p)))
  expect_gt(length(unique(vapply(sp$samples, function(s) max(s[[1]]), 1L))), 1L)
  emp <- as.numeric(table(factor(ks, levels = names(p)))) / length(ks)
  expect_lt(0.5 * sum(abs(emp - p)), 0.03)
})

test_that("nested weighted MCMC changes the top level and minimizers pad the hierarchy", {
  N <- 5L
  e <- rbind(c(1, 2), c(2, 3), c(1, 3), c(3, 4), c(4, 5), c(2, 5), c(1, 4))
  g <- graph_add_edges(graph_create(N, directed = FALSE), e)
  set.seed(3)
  w <- rexp(nrow(e))
  mk <- function(b0, b1) nested_weighted_block_state(
    g, weights = w, types = "real-exponential", levels = list(b0, b1), degree_corrected = FALSE
  )
  key <- function(b0, b1) paste(.canon(b0), .canon(b1[unique(b0)]), sep = "|")
  H <- list()
  for (b0 in .set_partitions(N)) for (b1 in .set_partitions(max(b0))) H[[length(H) + 1L]] <- list(b0, b1)
  S <- vapply(H, function(h) nested_weighted_description_length(mk(h[[1]], h[[2]])), numeric(1))
  p <- exp(-(S - min(S)))
  p <- p / sum(p)
  names(p) <- vapply(H, function(h) key(h[[1]], h[[2]]), "")
  set.seed(11)
  sp <- sample_posterior(mk(rep(1L, N), 1L), nsamples = 12000L, thin = 1L, burnin = 200L,
                         mcmc_args = list(d = 0.3))
  ks <- vapply(sp$samples, function(s) key(s[[1]], s[[2]]), "")
  emp <- as.numeric(table(factor(ks, levels = names(p)))) / length(ks)
  expect_lt(0.5 * sum(abs(emp - p)), 0.03)
  expect_gt(max(vapply(sp$samples, function(s) max(s[[2]]), 1L)), 1L)  # top level split

  # minimizers keep graph-tool's depth ceil(log2 N) + 1 (padding has zero DL)
  m <- minimize_nested_weighted_blockmodel_dl(g, weights = w, types = "real-exponential")
  expect_length(m$levels, ceiling(log2(N)) + 1L)
  lv <- m$levels
  expect_equal(nested_weighted_description_length(m),
               nested_weighted_description_length(nested_weighted_block_state(
                 g, weights = w, types = "real-exponential", levels = lv[seq_len(2L)])),
               tolerance = 1e-9)
})
