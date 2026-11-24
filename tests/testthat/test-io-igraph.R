test_that("vertex, graph, and general edge properties work and survive mutations", {
  g <- graph_create(3L, directed = TRUE)
  graph_set_graph_property(g, "title", "toy")
  graph_set_vertex_property(g, "name", c("a", "b", "c"))
  graph_set_vertex_property(g, "pos", rbind(c(0, 1), c(2, 3), c(4, 5)))

  # Add edges out of source-vertex order with both numeric and character properties
  graph_add_edges(
    g,
    rbind(c(3, 1), c(1, 2), c(2, 3)),
    properties = list(weight = c(30, 10, 20), label = c("e31", "e12", "e23"))
  )

  expect_equal(graph_graph_properties(g), "title")
  expect_equal(graph_graph_property(g, "title"), "toy")
  expect_equal(sort(graph_vertex_properties(g)), c("name", "pos"))
  expect_equal(graph_vertex_property(g, "name"), c("a", "b", "c"))
  expect_equal(sort(graph_edge_properties(g)), c("label", "weight"))

  # graph_edges() returns edges in vertex-major order: (1,2), (2,3), (3,1)
  expect_equal(unname(graph_edges(g)), rbind(c(1L, 2L), c(2L, 3L), c(3L, 1L)))
  expect_equal(graph_edge_property(g, "weight"), c(10, 20, 30))
  expect_equal(graph_edge_property(g, "label"), c("e12", "e23", "e31"))

  # Alias shares attribute environment mutations
  alias <- g
  graph_add_vertices(g, 1L)
  expect_equal(graph_vertex_property(alias, "name"), c("a", "b", "c", NA_character_))
  expect_equal(nrow(graph_vertex_property(alias, "pos")), 4L)
})

test_that("read_gt parses pinned graph-tool 3.9 directed and undirected .gt fixtures", {
  f_dir <- test_path("fixtures", "gt39_directed.gt")
  g <- read_gt(f_dir)

  expect_s3_class(g, "rgraphtool_graph")
  expect_true(g$directed)
  expect_equal(graph_num_vertices(g), 4)
  expect_equal(graph_num_edges(g), 5)
  expect_equal(
    unname(graph_edges(g)),
    rbind(c(1L, 2L), c(2L, 3L), c(2L, 3L), c(3L, 3L), c(3L, 1L))
  )

  # Graph properties
  expect_equal(graph_graph_property(g, "title"), "gt39 directed fixture")
  expect_identical(graph_graph_property(g, "flag"), TRUE)
  expect_identical(graph_graph_property(g, "count"), 42L)
  expect_equal(graph_graph_property(g, "big"), 1e12)
  expect_equal(graph_graph_property(g, "score"), 3.14159)
  expect_equal(graph_graph_property(g, "tags"), c("alpha", "beta"))
  expect_s3_class(graph_graph_property(g, "pyobj"), "gt_python_object")

  # Vertex properties
  expect_equal(graph_vertex_property(g, "name"), c("v1", "v2", "v3", "v4"))
  expect_identical(graph_vertex_property(g, "active"), c(TRUE, FALSE, TRUE, FALSE))
  expect_identical(graph_vertex_property(g, "s16"), c(-10L, 0L, 10L, 20L))
  expect_identical(graph_vertex_property(g, "i32"), c(100L, 200L, 300L, 400L))
  expect_identical(graph_vertex_property(g, "i64"), c(1L, 2L, 3L, 4L))
  expect_equal(graph_vertex_property(g, "val"), c(1.5, 2.5, 3.5, 4.5))
  expect_equal(graph_vertex_property(g, "pos"), list(c(0, 1), c(2, 3), c(4, 5), c(6, 7)))
  expect_equal(graph_vertex_property(g, "vbool"), list(c(TRUE, FALSE), logical(0), TRUE, c(FALSE, FALSE)))
  expect_equal(graph_vertex_property(g, "vint"), list(c(1L, 2L), 3L, integer(0), c(4L, 5L, 6L)))
  expect_equal(graph_vertex_property(g, "vstr"), list("a", c("b", "c"), character(0), "d"))

  # Edge properties (numeric ones usable directly in weighted_block_state)
  expect_equal(graph_edge_property(g, "weight"), c(1.1, 2.2, 3.3, 4.4, 5.5))
  expect_equal(graph_edge_property(g, "layer"), c(1, 2, 2, 3, 1))
  expect_equal(graph_edge_property(g, "label"), c("e01", "e12a", "e12b", "e22", "e20"))
  expect_identical(graph_edge_property(g, "is_loop"), c(FALSE, FALSE, FALSE, TRUE, FALSE))

  ws <- weighted_block_state(g, blocks = c(1L, 1L, 2L, 2L),
                             weights = "weight", types = "real-exponential")
  expect_true(is.finite(weighted_description_length(ws)))

  # Undirected gzip-compressed fixture
  f_undir <- test_path("fixtures", "gt39_undirected.gt.gz")
  gu <- read_gt(f_undir)
  expect_false(gu$directed)
  expect_equal(graph_num_vertices(gu), 4)
  expect_equal(graph_num_edges(gu), 4)
  expect_equal(
    unname(graph_edges(gu)),
    rbind(c(1L, 3L), c(2L, 2L), c(3L, 1L), c(4L, 2L))
  )
  expect_equal(graph_graph_property(gu, "name"), "gt39 undirected")
  expect_equal(graph_vertex_property(gu, "label"), c("a", "b", "c", "d"))
  expect_equal(graph_edge_property(gu, "weight"), c(10, 30, 20, 40))
})

test_that("write_gt and read_gt round-trip across compressions, edge cases, and N=255/256", {
  # Empty graph (N = 0)
  g0 <- graph_create(0L, directed = FALSE)
  f0 <- tempfile(fileext = ".gt")
  on.exit(unlink(f0), add = TRUE)
  write_gt(g0, f0)
  g0_in <- read_gt(f0)
  expect_false(g0_in$directed)
  expect_equal(graph_num_vertices(g0_in), 0)
  expect_equal(graph_num_edges(g0_in), 0)

  # Vertex index width boundary: N = 255 (uint8_t) vs N = 256 (uint16_t)
  for (n in c(255L, 256L)) {
    gn <- graph_create(n, directed = TRUE)
    graph_add_edges(gn, rbind(c(1L, n), c(n, 1L)))
    fn <- tempfile(fileext = ".gt")
    write_gt(gn, fn)
    gn_in <- read_gt(fn)
    unlink(fn)
    expect_equal(graph_num_vertices(gn_in), n)
    expect_equal(unname(graph_edges(gn_in)), rbind(c(1L, n), c(n, 1L)))
  }

  # Rich graph across uncompressed, gzip, bzip2, and xz
  g <- graph_create(4L, directed = TRUE)
  graph_add_edges(
    g,
    rbind(c(1, 2), c(2, 3), c(2, 3), c(3, 3)),
    properties = list(weight = c(0.5, 1.5, 2.5, 3.5),
                      tag = c("a", "b", "c", "d"),
                      ok = c(TRUE, FALSE, TRUE, FALSE))
  )
  graph_set_graph_property(g, "desc", "roundtrip")
  graph_set_graph_property(g, "nums", c(1.0, 2.0, 3.0))
  graph_set_vertex_property(g, "id", c(10L, 20L, 30L, 40L))
  graph_set_vertex_property(g, "coords", rbind(c(0, 1), c(2, 3), c(4, 5), c(6, 7)))

  for (ext in c(".gt", ".gt.gz", ".gt.bz2", ".gt.xz")) {
    fp <- tempfile(fileext = ext)
    write_gt(g, fp)
    g2 <- read_gt(fp)
    unlink(fp)
    expect_equal(graph_num_vertices(g2), 4)
    expect_equal(graph_num_edges(g2), 4)
    expect_equal(unname(graph_edges(g2)), unname(graph_edges(g)))
    expect_equal(graph_graph_property(g2, "desc"), "roundtrip")
    expect_equal(graph_graph_property(g2, "nums"), c(1.0, 2.0, 3.0))
    expect_identical(graph_vertex_property(g2, "id"), c(10L, 20L, 30L, 40L))
    expect_equal(graph_vertex_property(g2, "coords"),
                 list(c(0, 1), c(2, 3), c(4, 5), c(6, 7)))
    expect_equal(graph_edge_property(g2, "weight"), c(0.5, 1.5, 2.5, 3.5))
    expect_equal(graph_edge_property(g2, "tag"), c("a", "b", "c", "d"))
    expect_identical(graph_edge_property(g2, "ok"), c(TRUE, FALSE, TRUE, FALSE))
  }
})

test_that("read_gt decodes big-endian .gt files and rejects malformed files", {
  # Construct a minimal big-endian .gt file (N = 2, directed, 1 edge 0 -> 1, 1 int32 vertex prop)
  write_u64_be <- function(x) {
    # works for non-negative integers up to 2^31 - 1
    c(as.raw(rep(0L, 4L)), writeBin(as.integer(x), raw(), size = 4L, endian = "big"))
  }
  bytes <- c(
    as.raw(c(0xe2, 0x9b, 0xbe, 0x20, 0x67, 0x74)), # magic
    as.raw(0x01),                                   # version 1
    as.raw(0x01),                                   # big-endian
    write_u64_be(0L),                               # empty comment
    as.raw(0x01),                                   # directed
    write_u64_be(2L),                               # N = 2
    write_u64_be(1L), as.raw(0x01),                 # v=0: 1 neighbor -> 1
    write_u64_be(0L),                               # v=1: 0 neighbors
    write_u64_be(1L),                               # 1 property map
    as.raw(0x01),                                   # vertex prop
    write_u64_be(1L), charToRaw("x"),               # name "x"
    as.raw(0x02),                                   # int32_t
    writeBin(c(7L, 13L), raw(), size = 4L, endian = "big")
  )
  fb <- tempfile(fileext = ".gt")
  on.exit(unlink(fb), add = TRUE)
  writeBin(bytes, fb)

  gb <- read_gt(fb)
  expect_true(gb$directed)
  expect_equal(graph_num_vertices(gb), 2)
  expect_equal(unname(graph_edges(gb)), matrix(c(1L, 2L), nrow = 1L))
  expect_identical(graph_vertex_property(gb, "x"), c(7L, 13L))

  # Malformed magic header
  fbad <- tempfile(fileext = ".gt")
  on.exit(unlink(fbad), add = TRUE)
  writeBin(charToRaw("not a gt file"), fbad)
  expect_error(read_gt(fbad), "magic header")
})

test_that("as_igraph and from_igraph round-trip directed and undirected multigraphs with properties", {
  skip_if_not_installed("igraph")

  # Directed multigraph with parallel edges, self-loop, isolated vertex, and properties
  g <- graph_create(4L, directed = TRUE)
  graph_add_edges(
    g,
    rbind(c(1, 2), c(2, 3), c(2, 3), c(3, 3), c(3, 1)),
    properties = list(weight = c(1.0, 2.0, 3.0, 4.0, 5.0),
                      kind = c("a", "b", "c", "d", "e"))
  )
  graph_set_graph_property(g, "title", "igraph roundtrip")
  graph_set_vertex_property(g, "name", c("n1", "n2", "n3", "n4"))
  graph_set_vertex_property(g, "score", c(10, 20, 30, 40))

  ig <- as_igraph(g)
  expect_s3_class(ig, "igraph")
  expect_true(igraph::is_directed(ig))
  expect_equal(igraph::vcount(ig), 4L)
  expect_equal(igraph::ecount(ig), 5L)
  expect_equal(igraph::graph_attr(ig, "title"), "igraph roundtrip")
  expect_equal(igraph::vertex_attr(ig, "name"), c("n1", "n2", "n3", "n4"))
  expect_equal(igraph::vertex_attr(ig, "score"), c(10, 20, 30, 40))
  expect_equal(igraph::edge_attr(ig, "weight"), c(1.0, 2.0, 3.0, 4.0, 5.0))
  expect_equal(igraph::edge_attr(ig, "kind"), c("a", "b", "c", "d", "e"))

  g2 <- from_igraph(ig)
  expect_s3_class(g2, "rgraphtool_graph")
  expect_true(g2$directed)
  expect_equal(graph_num_vertices(g2), 4)
  expect_equal(graph_num_edges(g2), 5)
  expect_equal(unname(graph_edges(g2)), unname(graph_edges(g)))
  expect_equal(graph_graph_property(g2, "title"), "igraph roundtrip")
  expect_equal(graph_vertex_property(g2, "name"), c("n1", "n2", "n3", "n4"))
  expect_equal(graph_vertex_property(g2, "score"), c(10, 20, 30, 40))
  expect_equal(graph_edge_property(g2, "weight"), c(1.0, 2.0, 3.0, 4.0, 5.0))
  expect_equal(graph_edge_property(g2, "kind"), c("a", "b", "c", "d", "e"))

  # Native C++ SBM inference works directly on a graph imported from igraph
  ws <- weighted_block_state(g2, blocks = c(1L, 1L, 2L, 2L),
                             weights = "weight", types = "real-exponential")
  expect_true(is.finite(weighted_description_length(ws)))

  # Undirected igraph with out-of-order edges -> from_igraph -> write_gt -> read_gt
  igu <- igraph::make_graph(c(3, 1, 2, 2, 1, 3, 4, 2), n = 4, directed = FALSE)
  igraph::E(igu)$weight <- c(10, 20, 30, 40)
  igraph::E(igu)$tag <- c("e31", "e22", "e13", "e42")
  igraph::V(igu)$label <- c("v1", "v2", "v3", "v4")

  gu <- from_igraph(igu)
  expect_false(gu$directed)
  expect_equal(graph_num_vertices(gu), 4)
  expect_equal(graph_num_edges(gu), 4)

  # Check edge-to-property alignment matches igraph's edge table
  el_gt <- graph_edges(gu)
  w_gt <- graph_edge_property(gu, "weight")
  t_gt <- graph_edge_property(gu, "tag")
  # Edge (2,2) must have weight 20 and tag "e22"; edge (4,2) or (2,4) must have 40 and "e42"
  loop_idx <- which(el_gt[, 1] == 2L & el_gt[, 2] == 2L)
  expect_equal(w_gt[loop_idx], 20)
  expect_equal(t_gt[loop_idx], "e22")
  e42_idx <- which((el_gt[, 1] == 4L & el_gt[, 2] == 2L) |
                   (el_gt[, 1] == 2L & el_gt[, 2] == 4L))
  expect_equal(w_gt[e42_idx], 40)
  expect_equal(t_gt[e42_idx], "e42")

  # Empty igraph (n = 0)
  ig0 <- igraph::make_empty_graph(0L, directed = FALSE)
  g0 <- from_igraph(ig0)
  expect_equal(graph_num_vertices(g0), 0)
  expect_equal(graph_num_edges(g0), 0)
  ig0_back <- as_igraph(g0)
  expect_equal(igraph::vcount(ig0_back), 0L)
  expect_equal(igraph::ecount(ig0_back), 0L)
})

test_that("graph_load and graph_save support .gt, .gt.gz, and .graphml", {
  g <- graph_create(3L, directed = TRUE)
  graph_add_edges(g, rbind(c(1, 2), c(2, 3)), properties = list(weight = c(1.5, 2.5)))
  graph_set_vertex_property(g, "name", c("x", "y", "z"))

  f_gt <- tempfile(fileext = ".gt.gz")
  on.exit(unlink(f_gt), add = TRUE)
  graph_save(g, f_gt)
  g_gt <- graph_load(f_gt)
  expect_equal(graph_edge_property(g_gt, "weight"), c(1.5, 2.5))
  expect_equal(graph_vertex_property(g_gt, "name"), c("x", "y", "z"))

  skip_if_not_installed("igraph")
  f_gml <- tempfile(fileext = ".graphml")
  on.exit(unlink(f_gml), add = TRUE)
  graph_save(g, f_gml)
  g_gml <- graph_load(f_gml)
  expect_true(g_gml$directed)
  expect_equal(graph_num_vertices(g_gml), 3)
  expect_equal(graph_num_edges(g_gml), 2)
  expect_equal(graph_edge_property(g_gml, "weight"), c(1.5, 2.5))
  expect_equal(graph_vertex_property(g_gml, "name"), c("x", "y", "z"))
})
