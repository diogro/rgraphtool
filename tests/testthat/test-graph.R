test_that("directed multigraphs preserve loops and parallel edges", {
  g <- graph_create(3)
  graph_add_edges(g, rbind(c(1, 2), c(2, 2), c(2, 3), c(2, 3)))
  expect_equal(graph_num_vertices(g), 3)
  expect_equal(graph_num_edges(g), 4)
  expect_equal(unname(graph_edges(g)), rbind(c(1L, 2L), c(2L, 2L), c(2L, 3L), c(2L, 3L)))
  expect_equal(graph_degree(g, "in"), c(0, 2, 2))
  expect_equal(graph_degree(g, "out"), c(1, 3, 0))
  expect_equal(graph_degree(g), c(1, 5, 2))
})

test_that("undirected degrees, empty graphs, and aliases work", {
  g <- graph_create(directed = FALSE)
  expect_equal(dim(graph_edges(g)), c(0L, 2L))
  expect_length(graph_degree(g), 0)
  alias <- g
  graph_add_vertices(g, 2)
  graph_add_edges(g, rbind(c(2, 1), c(1, 1)))
  expect_equal(graph_num_edges(alias), 2)
  expect_equal(graph_degree(alias), c(3, 1))
  expect_error(graph_degree(g, "in"), "Undirected")
})

test_that("invalid inputs do not mutate a graph", {
  g <- graph_create(2)
  expect_error(graph_create(-1), "non-negative")
  expect_error(graph_create(directed = NA), "directed")
  expect_error(graph_add_vertices(g, 1.5), "non-negative")
  expect_error(graph_add_edges(g, rbind(c(1, 2), c(2, 3))), "valid 1-based")
  expect_error(graph_add_edges(g, matrix(c(1, NA_real_), ncol = 2)), "valid 1-based")
  expect_equal(graph_num_edges(g), 0)
  expect_error(graph_num_vertices(unserialize(serialize(g, NULL))), "cannot be restored")
})
