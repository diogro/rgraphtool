test_that("centrality, clustering, and modularity match pinned graph-tool 3.9 values", {
  # Two triangles connected by a bridge: (1,2,3) -- (4,5,6)
  g <- graph_create(6L, directed = FALSE)
  edges <- rbind(c(1, 2), c(2, 3), c(3, 1), c(3, 4), c(4, 5), c(5, 6), c(6, 4))
  w <- c(1.0, 2.0, 1.5, 0.5, 1.0, 2.5, 1.0)
  graph_add_edges(g, edges, properties = list(w = w))

  # PageRank (unweighted & weighted)
  expect_equal(
    as.numeric(pagerank(g)),
    c(0.1459854048150547, 0.1459854048150547, 0.20802919036989048,
      0.20802919036989048, 0.1459854048150547, 0.1459854048150547),
    tolerance = 1e-5
  )
  expect_equal(
    as.numeric(pagerank(g, weights = "w")),
    c(0.13673284247841072, 0.16001584613443584, 0.2082978127770684,
      0.13463701083097024, 0.18015824388955728, 0.18015824388955728),
    tolerance = 1e-5
  )

  # Brandes vertex & edge betweenness
  bw <- betweenness(g)
  expect_equal(bw$vertex, c(0, 0, 0.6, 0.6, 0, 0), tolerance = 1e-12)
  # Note: graph_edges(g) orders edges by (source, out_edges):
  # (1,2)[1.0], (1,3)[1.5], (2,3)[2.0], (3,4)[0.5], (4,5)[1.0], (4,6)[1.0], (5,6)[2.5]
  expect_equal(sort(bw$edge), sort(c(1/15, 4/15, 4/15, 9/15, 4/15, 1/15, 4/15)), tolerance = 1e-12)
  expect_equal(max(bw$edge), 0.6, tolerance = 1e-12)

  # Closeness & harmonic closeness
  expect_equal(closeness(g), c(0.5, 0.5, 5/7, 5/7, 0.5, 0.5), tolerance = 1e-12)
  expect_true(all(closeness(g, harmonic = TRUE) > 0))

  # Eigenvector, Katz, and HITS centrality
  ev <- eigenvector_centrality(g)
  expect_equal(ev$eigenvalue, 2.4142135623730754, tolerance = 1e-5)
  expect_equal(ev$vector, c(0.3535534, 0.3535534, 0.5, 0.5, 0.3535534, 0.3535534), tolerance = 1e-5)

  kc <- katz_centrality(g, alpha = 0.1)
  expect_equal(sum(kc^2), 1.0, tolerance = 1e-10)
  expect_gt(kc[3], kc[1])

  hc <- hits_centrality(g)
  expect_equal(hc$authority, hc$hub, tolerance = 1e-5)

  # Local & global clustering (with jackknife standard error matching graph-tool 3.9!)
  expect_equal(local_clustering(g), c(1, 1, 1/3, 1/3, 1, 1), tolerance = 1e-12)
  gc <- global_clustering(g)
  expect_equal(unname(gc["clustering"]), 0.6, tolerance = 1e-12)
  expect_equal(unname(gc["std_err"]), 0.18445509900106602, tolerance = 1e-12)

  # Modularity (gamma = 1, gamma = 2, and weighted)
  b <- c(1L, 1L, 1L, 2L, 2L, 2L)
  expect_equal(modularity(g, b), 5 / 14, tolerance = 1e-12)
  expect_equal(modularity(g, b, gamma = 2.0), -2 / 14, tolerance = 1e-12)
  expect_equal(modularity(g, b, weights = "w"), 8.5 / 19.0, tolerance = 1e-12)
})

test_that("shortest paths, diameter, components, biconnectivity, DAG, and assortativity work", {
  # Barbell graph: bridge is edge (3, 4), articulation points are 3 and 4
  g <- graph_create(6L, directed = FALSE)
  graph_add_edges(g, rbind(c(1, 2), c(2, 3), c(1, 3), c(3, 4), c(4, 5), c(5, 6), c(4, 6)))

  dmat <- shortest_distance(g)
  expect_equal(dim(dmat), c(6L, 6L))
  expect_equal(dmat[1, 6], 3)
  expect_equal(graph_diameter(g), 3)
  expect_equal(pseudo_diameter(g)$diameter, 3)

  sp <- shortest_path(g, 1L, 6L)
  expect_equal(sp$distance, 3)
  expect_equal(sp$vertices, c(1L, 3L, 4L, 6L))
  expect_length(sp$edges, 3L)

  bic <- label_biconnected_components(g)
  expect_equal(bic$num_components, 3L)
  expect_equal(sort(bic$articulation_points), c(3L, 4L))
  expect_length(bic$bridges, 1L)

  # Directed SCCs, attractors, DAG, and Bellman-Ford negative weights
  gd <- graph_create(5L, directed = TRUE)
  # SCC {1, 2} -> vertex 3 -> sink SCC {4, 5}
  graph_add_edges(gd, rbind(c(1, 2), c(2, 1), c(2, 3), c(3, 4), c(4, 5), c(5, 4)))
  lc <- label_components(gd)
  expect_equal(lc$num_components, 3L)
  expect_equal(sum(lc$attractors), 1L)
  expect_false(is_dag(gd))
  expect_error(topological_sort(gd), "cycle")

  dag <- graph_create(4L, directed = TRUE)
  graph_add_edges(dag, rbind(c(1, 2), c(1, 3), c(2, 4), c(3, 4)),
                  properties = list(w = c(4.0, 1.0, -2.0, 2.0)))
  expect_true(is_dag(dag))
  ord <- topological_sort(dag)
  expect_equal(ord[1L], 1L)
  expect_equal(ord[4L], 4L)

  # Bellman-Ford picks 1 -> 2 -> 4 (cost 4 + (-2) = 2) over 1 -> 3 -> 4 (cost 3)
  sp_neg <- shortest_path(dag, 1L, 4L, weights = "w")
  expect_equal(sp_neg$distance, 2.0)
  expect_equal(sp_neg$vertices, c(1L, 2L, 4L))

  # Assortativity (categorical & scalar)
  expect_gt(assortativity(g, c("A", "A", "A", "B", "B", "B"), categorical = TRUE), 0.7)
  expect_true(is.finite(assortativity(g, "total")))
})

test_that("random graph generators and modularity overfitting warning behave as expected", {
  set.seed(42)

  g_gnp <- random_graph_erdos_renyi(100L, p = 0.06, directed = FALSE)
  expect_equal(graph_num_vertices(g_gnp), 100)

  g_gnm <- random_graph_erdos_renyi(50L, m = 120L, directed = TRUE)
  expect_equal(graph_num_edges(g_gnm), 120)

  g_ba <- random_graph_barabasi_albert(80L, m = 2L)
  expect_equal(graph_num_vertices(g_ba), 80)
  expect_gt(graph_num_edges(g_ba), 100)

  g_ws <- random_graph_watts_strogatz(60L, k = 4L, p = 0.1)
  expect_equal(graph_num_vertices(g_ws), 60)
  expect_equal(graph_num_edges(g_ws), 120)

  g_cfg <- random_graph_configuration(rep(4L, 30L))
  expect_equal(graph_degree(g_cfg), rep(4, 30))

  P <- rbind(c(0.6, 0.02), c(0.02, 0.6))
  g_sbm <- random_graph_sbm(rep(1:2, each = 25L), P)
  expect_equal(graph_num_vertices(g_sbm), 50)

  # Peixoto demonstration: maximize_modularity MUST emit a loud warning by default
  # and overfits a pure Erdos-Renyi random graph into multiple spurious communities,
  # whereas MDL inference (minimize_blockmodel_dl) rightly infers B = 1!
  expect_warning(mod_er <- maximize_modularity(g_gnp), "terrible idea")
  expect_s3_class(mod_er, "rgraphtool_modularity_state")
  expect_gt(mod_er$B, 3L)
  expect_gt(mod_er$modularity, 0.25)

  sbm_er <- minimize_blockmodel_dl(g_gnp)
  expect_equal(length(block_sizes(sbm_er)), 1L)

  # On a planted SBM, warn = FALSE suppresses warning and recovers planted structure
  expect_silent(mod_sbm <- maximize_modularity(g_sbm, warn = FALSE))
  expect_equal(mod_sbm$B, 2L)
  expect_equal(partition_overlap(block_membership(mod_sbm), rep(1:2, each = 25L)), 1.0)
})
