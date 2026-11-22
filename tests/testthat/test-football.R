make_football_archive <- function(bad_position = FALSE, bad_edge = FALSE) {
  path <- tempfile()
  dir.create(path)
  nodes <- data.frame(index = 114:0, label = paste0("team", 114:0),
                      conference = rep(0L, 115),
                      position = sprintf("array([%d, %d])", 114:0, -(114:0)))
  if (bad_position) nodes$position[1] <- "array([oops, 0])"
  edges <- as.data.frame(t(combn(0:114, 2)))[seq_len(613), ]
  if (bad_edge) edges[1, 1] <- 115
  utils::write.table(nodes, file.path(path, "nodes.csv"), sep = ",",
                     row.names = FALSE, col.names = FALSE)
  utils::write.table(edges, file.path(path, "edges.csv"), sep = ",",
                     row.names = FALSE, col.names = FALSE)
  archive <- file.path(path, "football.zip")
  utils::zip(archive, file.path(path, c("nodes.csv", "edges.csv")), flags = "-q -j")
  archive
}

test_that("football archives load in vertex order with one-based edges", {
  skip_if(Sys.which("zip") == "")
  archive <- make_football_archive()
  on.exit(unlink(dirname(archive), recursive = TRUE))
  football <- football_network(archive)
  expect_equal(graph_num_vertices(football$graph), 115)
  expect_equal(graph_num_edges(football$graph), 613)
  expect_equal(football$nodes$index, 1:115)
  expect_equal(football$nodes$label[1], "team0")
  expect_equal(unname(football$positions[2, ]), c(1, -1))
  expect_equal(unname(graph_edges(football$graph)[1, ]), c(1L, 2L))
})

test_that("invalid football inputs fail explicitly", {
  skip_if(Sys.which("zip") == "")
  expect_error(football_network("missing-football.zip"), "path")
  bad_position <- make_football_archive(bad_position = TRUE)
  bad_edge <- make_football_archive(bad_edge = TRUE)
  on.exit(unlink(c(dirname(bad_position), dirname(bad_edge)), recursive = TRUE))
  expect_error(football_network(bad_position), "non-numeric node positions")
  expect_error(football_network(bad_edge), "valid 1-based")
})
