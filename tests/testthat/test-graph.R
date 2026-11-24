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

test_that("print.rgraphtool_graph formats small, large, and deserialized graphs", {
  g <- graph_create(4L, directed = TRUE)
  graph_add_edges(g, rbind(c(1, 2), c(2, 3), c(2, 3), c(3, 3), c(3, 1)),
                  properties = list(weight = c(1, 2, 3, 4, 5)))
  graph_set_vertex_property(g, "name", c("a", "b", "c", "d"))
  graph_set_graph_property(g, "title", "demo")

  out <- capture.output(print(g))
  expect_true(any(grepl("rgraphtool graph \\(directed, 4 vertices, 5 edges\\)", out)))
  expect_true(any(grepl("graph props: title \\(chr\\)", out)))
  expect_true(any(grepl("vertex props: name \\(chr\\)", out)))
  expect_true(any(grepl("edge props: weight \\(dbl\\)", out)))
  expect_true(any(grepl("adjacency:", out)))
  expect_true(any(grepl("a->b", out)))

  # Large graph skips ASCII adjacency by default and truncates edges
  gl <- graph_create(20L, directed = FALSE)
  graph_add_edges(gl, cbind(1:19, 2:20))
  out_l <- capture.output(print(gl, max_edges = 5L))
  expect_false(any(grepl("adjacency:", out_l)))
  expect_true(any(grepl("\\+14 more", out_l)))

  # Deserialized pointer prints gracefully instead of erroring
  dead <- unserialize(serialize(g, NULL))
  out_dead <- capture.output(print(dead))
  expect_true(any(grepl("invalid or deserialized", out_dead)))
})
