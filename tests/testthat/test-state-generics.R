test_that("description_length() and state_blocks() dispatch on every state family", {
  set.seed(3)
  g <- graph_create(12, directed = FALSE)
  graph_add_edges(g, rbind(t(utils::combn(1:6, 2)), t(utils::combn(7:12, 2)), c(1, 7)),
                  properties = list(w = stats::rnorm(31), layer = rep(1:2, length.out = 31)))
  b <- rep(1:2, each = 6)
  lv <- list(b, c(1L, 1L))
  states <- list(
    flat = list(block_state(g, b), block_description_length),
    nested = list(nested_block_state(g, lv), nested_description_length),
    weighted = list(weighted_block_state(g, "w", "real-normal", blocks = b),
                    weighted_description_length),
    nested_weighted = list(nested_weighted_block_state(g, "w", "real-normal", levels = lv),
                           nested_weighted_description_length),
    layered = list(layered_block_state(g, "layer", blocks = b), layered_description_length),
    nested_layered = list(nested_layered_block_state(g, "layer", levels = lv),
                          nested_layered_description_length),
    pp = list(pp_block_state(g, b), pp_description_length)
  )
  for (nm in names(states)) {
    st <- states[[nm]][[1]]
    expect_equal(description_length(st), states[[nm]][[2]](st), info = nm)
    expect_equal(as.integer(state_blocks(st)), b, info = nm)
    expect_identical(is_nested_state(st), grepl("nested", nm), info = nm)
  }
  expect_equal(state_blocks(states$nested[[1]], level = 1), rep(1L, 12))
  expect_equal(description_length(states$flat[[1]], degree_prior = "uniform"),
               block_description_length(states$flat[[1]], degree_prior = "uniform"))
  expect_error(state_blocks(states$flat[[1]], level = 1), "only level 0")
  expect_error(description_length(list()), "no description length")
})

test_that("is_block_state() recognises states and state_blocks() projects dynamics hierarchies", {
  g <- graph_create(4, directed = FALSE)
  graph_add_edges(g, rbind(c(1, 2), c(3, 4)))
  expect_true(is_block_state(block_state(g, c(1, 1, 2, 2))))
  expect_false(is_block_state(g))
  expect_false(is_block_state(1:4))
  fake <- structure(list(bs = list(c(1L, 1L, 2L, 3L), c(1L, 2L, 2L), c(1L, 1L))),
                    class = "rgraphtool_dynamics_block_state")
  expect_true(is_block_state(fake))
  expect_equal(state_blocks(fake), c(1L, 1L, 2L, 3L))
  expect_equal(state_blocks(fake, 1), c(1L, 1L, 2L, 2L))
  expect_error(state_blocks(fake, 3), "between 0 and 2")
})
