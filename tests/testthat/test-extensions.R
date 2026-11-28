test_that("kcore_decomposition matches graph-tool 3.9 across checked and unchecked modes", {
  # Empty and single-vertex graphs
  expect_equal(kcore_decomposition(graph_create(0L, directed = FALSE)), integer(0))
  expect_equal(kcore_decomposition(graph_create(1L, directed = FALSE)), 0L)

  # Reference fixture verified against pinned graph-tool 3.9:
  # 4-clique (1..4) + parallel edge (4,5) + edge (5,6) + self-loop (6,6) + isolated 7,8
  g <- graph_create(8L, directed = FALSE)
  graph_add_edges(g, rbind(
    c(1L, 2L), c(2L, 3L), c(3L, 1L),
    c(1L, 4L), c(2L, 4L), c(3L, 4L),
    c(4L, 5L), c(5L, 6L), c(6L, 6L), c(4L, 5L)
  ))
  kc_u <- kcore_decomposition(g, checked = FALSE)
  kc_c <- kcore_decomposition(g, checked = TRUE)
  expect_equal(kc_u, c(3L, 3L, 3L, 3L, 3L, 3L, 0L, 0L))
  expect_equal(kc_c, kc_u)

  # Directed reference fixture verified against pinned graph-tool 3.9
  gd <- graph_create(6L, directed = TRUE)
  graph_add_edges(gd, rbind(
    c(1L, 2L), c(2L, 3L), c(3L, 1L),
    c(1L, 3L), c(2L, 1L), c(3L, 2L),
    c(3L, 4L), c(4L, 5L), c(5L, 4L)
  ))
  expect_equal(kcore_decomposition(gd, deg = "total", checked = FALSE), c(4L, 4L, 4L, 2L, 2L, 0L))
  expect_equal(kcore_decomposition(gd, deg = "total", checked = TRUE),  c(4L, 4L, 4L, 2L, 2L, 0L))

  # Directed out-core and in-core on asymmetric DAG + cycle
  gd2 <- graph_create(5L, directed = TRUE)
  graph_add_edges(gd2, rbind(
    c(1L, 2L), c(2L, 3L), c(3L, 1L), # 1-cycle
    c(4L, 1L), c(4L, 2L),            # 4 has out-deg 2, in-deg 0
    c(3L, 5L)                        # 5 has out-deg 0, in-deg 1
  ))
  expect_equal(kcore_decomposition(gd2, deg = "out", checked = FALSE),
               kcore_decomposition(gd2, deg = "out", checked = TRUE))
  expect_equal(kcore_decomposition(gd2, deg = "in", checked = FALSE),
               kcore_decomposition(gd2, deg = "in", checked = TRUE))
  expect_equal(kcore_decomposition(gd2, deg = "out")[5L], 0L)
  expect_equal(kcore_decomposition(gd2, deg = "in")[4L], 0L)

  expect_error(kcore_decomposition(g, deg = "out"), "Undirected")
})

test_that("vertex_strength reads attached edge properties identically in checked/unchecked modes", {
  g <- graph_create(5L, directed = TRUE)
  edges <- rbind(c(1L, 2L), c(2L, 3L), c(1L, 3L), c(3L, 3L), c(1L, 2L))
  w <- c(1.5, 2.0, 0.5, 4.0, 3.5)
  graph_add_edges(g, edges, properties = list(weight = w))

  s_out_u <- vertex_strength(g, "weight", mode = "out", checked = FALSE)
  s_out_c <- vertex_strength(g, "weight", mode = "out", checked = TRUE)
  expect_equal(s_out_u, c(5.5, 2.0, 4.0, 0.0, 0.0))
  expect_equal(s_out_c, s_out_u)

  s_in_u <- vertex_strength(g, "weight", mode = "in", checked = FALSE)
  expect_equal(s_in_u, c(0.0, 5.0, 6.5, 0.0, 0.0))

  s_tot_u <- vertex_strength(g, "weight", mode = "total", checked = FALSE)
  expect_equal(s_tot_u, s_out_u + s_in_u)

  expect_error(vertex_strength(g, "missing_prop"), "does not exist")
})

test_that("checked vs unchecked property map bounds checks work as documented", {
  g <- graph_create(4L, directed = FALSE)
  vprop <- c(10, 20, 30, 40)

  # Valid checked and unchecked reads
  expect_equal(rgraphtool:::cpp_extension_check_bounds(g, vprop, 0L, FALSE), 10)
  expect_equal(rgraphtool:::cpp_extension_check_bounds(g, vprop, 3L, TRUE), 40)

  # Out-of-bounds checked read throws std::out_of_range
  expect_error(rgraphtool:::cpp_extension_check_bounds(g, vprop, 4L, FALSE), "out of bounds")

  # Creating an unchecked view over a too-short vector throws up-front
  short_vec <- c(1, 2)
  expect_error(rgraphtool:::cpp_extension_check_bounds(g, short_vec, 0L, TRUE), "smaller than num_vertices")

  # checked_vertex rejects out-of-range vertex index
  expect_error(rgraphtool:::cpp_extension_check_bounds(g, vprop, 9L, TRUE), "Vertex index out of bounds")
})

test_that("bundled C++ extension examples and inline user C++ code compile and run via compile_cpp_extension", {
  inc_dir <- rgraphtool_include_dir()
  expect_true(file.exists(file.path(inc_dir, "rgraphtool.h")))

  kcore_cpp <- rgraphtool_extension_example("kcore")
  wdeg_cpp <- rgraphtool_extension_example("weighted_degree")
  expect_true(file.exists(kcore_cpp))
  expect_true(file.exists(wdeg_cpp))

  # Compile and test standalone k-core extension
  ext_env <- new.env(parent = globalenv())
  compile_cpp_extension(file = kcore_cpp, env = ext_env)
  expect_true(exists("ext_kcore_decomposition", envir = ext_env))

  g <- graph_create(8L, directed = FALSE)
  graph_add_edges(g, rbind(
    c(1L, 2L), c(2L, 3L), c(3L, 1L),
    c(1L, 4L), c(2L, 4L), c(3L, 4L),
    c(4L, 5L), c(5L, 6L), c(6L, 6L), c(4L, 5L)
  ), properties = list(weight = seq(0.5, 5.0, by = 0.5)))

  expect_equal(ext_env$ext_kcore_decomposition(g, checked = FALSE),
               kcore_decomposition(g, checked = FALSE))
  expect_equal(ext_env$ext_kcore_decomposition(g, checked = TRUE),
               kcore_decomposition(g, checked = TRUE))

  # Compile and test standalone weighted-strength extension
  compile_cpp_extension(file = wdeg_cpp, env = ext_env)
  expect_true(exists("ext_weighted_strength", envir = ext_env))
  expect_equal(ext_env$ext_weighted_strength(g, "weight", checked = FALSE),
               vertex_strength(g, "weight", checked = FALSE))
  expect_equal(ext_env$ext_weighted_strength(g, "weight", checked = TRUE),
               vertex_strength(g, "weight", checked = TRUE))

  # Compile inline user C++ extension exercising edges_range, out_edges_range,
  # in_edges_range, make_edge_property, and checked EdgePropertyMap error
  inline_code <- '
#include <rgraphtool.h>

// [[Rcpp::export]]
Rcpp::List ext_inspect_graph(SEXP graph) {
    const rgraphtool::Graph& g = rgraphtool::get_graph(graph);
    auto ep = rgraphtool::make_edge_property<int>(g, 0);
    auto ep_u = ep.get_unchecked();
    int total_e = 0;
    for (const auto& e : rgraphtool::edges_range(g)) {
        ep_u[e] = static_cast<int>(rgraphtool::source(e, g) + rgraphtool::target(e, g));
        total_e++;
    }
    int sum_ep = 0;
    for (const auto& e : rgraphtool::edges_range(g)) {
        sum_ep += ep[e];
    }
    bool caught_oob = false;
    try {
        (void) ep[rgraphtool::num_edges(g) + 5];
    } catch (const std::out_of_range&) {
        caught_oob = true;
    }
    return Rcpp::List::create(
        Rcpp::Named("n") = static_cast<int>(rgraphtool::num_vertices(g)),
        Rcpp::Named("m") = total_e,
        Rcpp::Named("sum_ep") = sum_ep,
        Rcpp::Named("caught_oob") = caught_oob
    );
}
'
  compile_cpp_extension(code = inline_code, env = ext_env)
  info <- ext_env$ext_inspect_graph(g)
  expect_equal(info$n, 8L)
  expect_equal(info$m, 10L)
  expect_true(info$caught_oob)
  expect_gt(info$sum_ep, 0L)
})
