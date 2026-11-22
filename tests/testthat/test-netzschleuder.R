make_archive <- function(directed = TRUE, header = TRUE) {
  path <- tempfile()
  dir.create(path)
  write_csv <- function(lines, name) writeLines(lines, file.path(path, name))
  write_csv(c(if (header) "# index, label, _pos",
              '1,b,"array([1.5, -2])"', '0,a,"array([0, 1])"', '2,c,"array([3, 4])"'),
            "nodes.csv")
  write_csv(c(if (header) "# source, target, value", "0,1,5", "1,2,7", "2,2,1"),
            "edges.csv")
  write_csv(c("# prop_name, value", "name,toy", sprintf("directed,%d.0", directed)),
            "gprops.csv")
  archive <- file.path(path, "toy.zip")
  utils::zip(archive, file.path(path, c("nodes.csv", "edges.csv", "gprops.csv")),
             flags = "-q -j")
  archive
}

test_that("Netzschleuder archives load with attributes and positions", {
  skip_if(Sys.which("zip") == "")
  archive <- make_archive()
  on.exit(unlink(dirname(archive), recursive = TRUE))
  net <- netzschleuder_network("toy", file = archive)
  expect_equal(graph_num_vertices(net$graph), 3)
  expect_equal(graph_num_edges(net$graph), 3)
  expect_equal(graph_degree(net$graph, "out"), c(1, 1, 1))
  expect_equal(net$nodes$label, c("a", "b", "c"))
  expect_equal(net$nodes$index, 1:3)
  expect_equal(unname(net$positions[2, ]), c(1.5, -2))
  expect_equal(net$edges$value, c(5, 7, 1))
  expect_equal(net$edges$source, c(1, 2, 3))
  expect_equal(net$properties$name, "toy")
  undirected <- netzschleuder_network("toy", file = archive, directed = FALSE)
  expect_error(graph_degree(undirected$graph, "out"), "Undirected")
})

test_that("Netzschleuder inputs are validated", {
  skip_if(Sys.which("zip") == "")
  expect_error(netzschleuder_network("bad name"), "identifiers")
  expect_error(netzschleuder_network("toy", file = "missing.zip"), "path")
  archive <- make_archive(header = FALSE)
  on.exit(unlink(dirname(archive), recursive = TRUE))
  expect_error(netzschleuder_network("toy", file = archive), "header")
})
