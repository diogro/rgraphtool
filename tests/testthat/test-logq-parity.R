# Degree DLs with blocks holding more than 2000 edge ends, pinned to
# graph-tool 3.9 (src/logq.h). Model-specific log q(n, k) copies used to cap
# the exact table at 2000 and fall back to an incorrect approximation
# (hundreds to thousands of nats off).
.logq_fixture <- function() {
  fx <- function(x) test_path("fixtures", sprintf("logq_parity_%s.csv", x))
  e <- as.matrix(utils::read.csv(fx("e"), header = FALSE)) + 1L
  w <- utils::read.csv(fx("w"), header = FALSE)[[1L]]
  b <- utils::read.csv(fx("b"), header = FALSE)[[1L]] + 1L
  g <- graph_create(length(b), directed = FALSE)
  graph_add_edges(g, e, properties = list(w = w))
  list(g = g, b = b, w = graph_edge_property(g, "w"))
}

test_that("weighted DC degree DL matches graph-tool beyond the old 2000 cap", {
  fx <- .logq_fixture()
  ref <- c(distributed = 14524.5470885429, uniform = 14706.9295788729)
  for (kind in names(ref)) {
    st <- weighted_block_state(fx$g, weights = fx$w, types = "real-exponential",
                               blocks = fx$b, degree_corrected = TRUE, degree_prior = kind)
    expect_equal(weighted_description_length(st), ref[[kind]], tolerance = 1e-8)
  }
})

test_that("layered and planted-partition DC degree DLs match graph-tool", {
  fx <- .logq_fixture()
  E <- graph_num_edges(fx$g)
  ref_l <- c(distributed = 11937.905908, uniform = 12120.288398)
  ref_pp <- c(distributed = 11937.946714, uniform = 12120.329204)
  for (kind in names(ref_l)) {
    ls <- layered_block_state(fx$g, layers = rep(1L, E), blocks = fx$b, independent = TRUE,
                              degree_corrected = TRUE, degree_prior = kind)
    expect_equal(layered_description_length(ls), ref_l[[kind]], tolerance = 1e-6)
    pp <- pp_block_state(fx$g, blocks = fx$b, degree_prior = kind)
    expect_equal(pp_description_length(pp), ref_pp[[kind]], tolerance = 1e-6)
  }
})
