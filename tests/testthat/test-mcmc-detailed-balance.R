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
