test_that("sfdp_layout handles edge cases, determinism, pinning, groups, and quadtree", {
  # Empty graph
  g0 <- graph_create(0L, directed = FALSE)
  p0 <- sfdp_layout(g0)
  expect_equal(dim(p0), c(0L, 2L))
  expect_equal(colnames(p0), c("x", "y"))

  # Single vertex
  g1 <- graph_create(1L, directed = FALSE)
  p1 <- sfdp_layout(g1)
  expect_equal(dim(p1), c(1L, 2L))
  expect_true(all(is.finite(p1)))

  # Determinism under set.seed() and self-loops/parallel edges/disconnected nodes
  g <- graph_create(12L, directed = FALSE)
  edges <- rbind(
    cbind(1:5, c(2:5, 1)),
    cbind(6:10, c(7:10, 6)),
    c(1L, 1L),   # self-loop
    c(1L, 2L)    # parallel edge
  )
  graph_add_edges(g, edges) # vertices 11 and 12 are isolated
  set.seed(42)
  pos_a <- sfdp_layout(g, max_iter = 40L)
  set.seed(42)
  pos_b <- sfdp_layout(g, max_iter = 40L)
  expect_equal(pos_a, pos_b)
  expect_equal(dim(pos_a), c(12L, 2L))
  expect_true(all(is.finite(pos_a)))

  # Warm-start + pinned vertices remain fixed
  pin_mask <- rep(FALSE, 12L)
  pin_mask[c(1L, 6L)] <- TRUE
  set.seed(7)
  pos_warm <- sfdp_layout(g, pos = pos_a, init_step = 0.05, max_iter = 5L, pin = pin_mask)
  expect_equal(pos_warm[c(1L, 6L), ], pos_a[c(1L, 6L), ])
  expect_gt(max(abs(pos_warm[-c(1L, 6L), ] - pos_a[-c(1L, 6L), ])), 0)

  # Group attraction pulls intra-group vertices closer
  g_grp <- graph_create(20L, directed = FALSE)
  graph_add_edges(g_grp, cbind(1:19, 2:20))
  grp <- rep(1:2, each = 10L)
  set.seed(99)
  pos_g <- sfdp_layout(g_grp, groups = grp, mu = 2.0, gamma = -0.5, max_iter = 80L)
  d_mat <- as.matrix(stats::dist(pos_g))
  intra <- mean(d_mat[1:10, 1:10][upper.tri(d_mat[1:10, 1:10])])
  inter <- mean(d_mat[1:10, 11:20])
  expect_lt(intra, inter)

  # QuadTree path (N > 128)
  g_big <- graph_create(160L, directed = FALSE)
  ring_e <- cbind(seq_len(160L), c(seq_len(159L) + 1L, 1L))
  graph_add_edges(g_big, ring_e)
  set.seed(11)
  pos_big <- sfdp_layout(g_big, max_iter = 25L)
  expect_equal(dim(pos_big), c(160L, 2L))
  expect_true(all(is.finite(pos_big)))

  # Invalid inputs
  expect_error(sfdp_layout(g, K = -1), "positive")
  expect_error(sfdp_layout(g, pos = matrix(1:4, 2, 2)), "12 x 2")
  expect_error(sfdp_layout(g, groups = 1:3), "groups")
  expect_error(sfdp_layout(g, pin = TRUE), "pin")
})

test_that("graph_bfs computes exact distances, parents, and tree edges", {
  # Directed path 1 -> 2 -> 3 -> 4 plus branch 2 -> 5, isolated 6, parallel edge, self-loop
  g <- graph_create(6L, directed = TRUE)
  e <- rbind(
    c(1L, 2L),
    c(2L, 3L),
    c(3L, 4L),
    c(2L, 5L),
    c(2L, 3L), # parallel edge
    c(1L, 1L)  # self-loop
  )
  graph_add_edges(g, e)

  b_dir <- graph_bfs(g, source = 2L, directed = TRUE)
  expect_s3_class(b_dir, "rgraphtool_bfs")
  expect_equal(b_dir$distance, c(Inf, 0, 1, 2, 1, Inf))
  expect_true(is.na(b_dir$parent[2L]))
  expect_true(is.na(b_dir$parent[1L]))
  expect_true(is.na(b_dir$parent[6L]))
  expect_equal(b_dir$parent[c(3L, 4L, 5L)], c(2L, 3L, 2L))
  expect_equal(sum(b_dir$tree_edges), 3L) # parallel edge counted once

  # Undirected traversal from vertex 2 reaches vertex 1 too
  b_undir <- graph_bfs(g, source = 2L, directed = FALSE)
  expect_equal(b_undir$distance, c(1, 0, 1, 2, 1, Inf))
  expect_equal(sum(b_undir$tree_edges), 4L)

  # max_dist truncation
  b_lim <- graph_bfs(g, source = 1L, max_dist = 2L, directed = TRUE)
  expect_equal(b_lim$distance, c(0, 1, 2, Inf, 2, Inf))
  expect_equal(sum(b_lim$tree_edges), 3L)

  b_zero <- graph_bfs(g, source = 1L, max_dist = 0L, directed = TRUE)
  expect_equal(b_zero$distance, c(0, Inf, Inf, Inf, Inf, Inf))
  expect_equal(sum(b_zero$tree_edges), 0L)

  expect_output(print(b_dir), "rgraphtool BFS traversal")

  # Invalid inputs
  g0 <- graph_create(0L, directed = FALSE)
  expect_error(graph_bfs(g0, 1L), "no vertices")
  expect_error(graph_bfs(g, 0L), "source")
  expect_error(graph_bfs(g, 7L), "source")
  expect_error(graph_bfs(g, 1L, max_dist = -1L), "max_dist")
})

test_that("plot_graph renders undirected and directed graphs cleanly offscreen", {
  g <- graph_create(8L, directed = TRUE)
  graph_add_edges(g, cbind(c(1, 2, 3, 4, 5, 6, 7, 1), c(2, 3, 4, 5, 6, 7, 8, 5)))
  set.seed(123)
  pos <- sfdp_layout(g, max_iter = 30L)

  png_file <- tempfile("rgt_plot_", fileext = ".png")
  on.exit(unlink(png_file), add = TRUE)
  grDevices::png(png_file, width = 320, height = 240)
  out_pos <- tryCatch(
    plot_graph(
      g, pos = pos,
      vertex_color = rep(c("#3B82F6", "#EF4444"), 4L),
      vertex_halo = c(TRUE, FALSE, TRUE, rep(FALSE, 5L)),
      vertex_visible = c(rep(TRUE, 7L), FALSE),
      vertex_label = seq_len(8L),
      main = "Test graph"
    ),
    finally = grDevices::dev.off()
  )
  expect_equal(out_pos, pos)
  expect_true(file.exists(png_file))
  expect_gt(file.info(png_file)$size, 500)
})

test_that("animate_sirs obeys exact SIRS mechanics and exports HTML and PNG frames", {
  # Small ring graph
  n <- 16L
  g <- graph_create(n, directed = FALSE)
  graph_add_edges(g, cbind(seq_len(n), c(seq_len(n - 1L) + 1L, 1L)))
  set.seed(5)
  pos <- sfdp_layout(g, max_iter = 30L)

  # Deterministic absorbing boundary checks
  s_zero <- animate_sirs(g, pos = pos, n_steps = 5L, x = 0, r = 0, s = 0,
                         init_state = rep(0L, n), mode = "none")
  expect_true(all(s_zero$states == 0L))
  expect_equal(as.integer(rowSums(s_zero$counts)), rep(n, 6L))

  s_rec <- animate_sirs(g, pos = pos, n_steps = 1L, x = 0, r = 1, s = 0,
                        init_state = rep(1L, n), mode = "none")
  expect_true(all(s_rec$states[, 2L] == 2L))

  s_res <- animate_sirs(g, pos = pos, n_steps = 1L, x = 0, r = 0, s = 1,
                        init_state = rep(2L, n), mode = "none")
  expect_true(all(s_res$states[, 2L] == 0L))

  # Fear state marks susceptible neighbors of infected vertex 1 (vertices 2 and 16 on ring)
  init_one <- integer(n)
  init_one[1L] <- 1L
  s_fear <- animate_sirs(g, pos = pos, n_steps = 0L, x = 0, r = 0, s = 0,
                         init_state = init_one, fear_state = TRUE, mode = "none")
  expect_equal(which(s_fear$fear[, 1L]), c(2L, 16L))

  # HTML export & reproducibility
  html_file <- tempfile("rgt_sirs_", fileext = ".html")
  on.exit(unlink(html_file), add = TRUE)
  set.seed(77)
  anim_html <- animate_sirs(g, pos = pos, n_steps = 20L, x = 0.05, r = 0.2, s = 0.05,
                            init_state = init_one, fear_state = TRUE,
                            mode = "html", output = html_file)
  expect_s3_class(anim_html, "rgraphtool_sirs_animation")
  expect_true(file.exists(html_file))
  txt <- paste(readLines(html_file, warn = FALSE), collapse = "\n")
  expect_match(txt, "<canvas id=\"cv\"")
  expect_match(txt, "requestAnimationFrame")
  expect_output(print(anim_html), "rgraphtool SIRS epidemic animation")

  # Reproducibility with same seed
  set.seed(77)
  anim_rep <- animate_sirs(g, pos = pos, n_steps = 20L, x = 0.05, r = 0.2, s = 0.05,
                           init_state = init_one, fear_state = TRUE, mode = "none")
  expect_equal(anim_html$states, anim_rep$states)
  expect_equal(anim_html$newly_infected, anim_rep$newly_infected)

  # PNG frames export
  frames_dir <- tempfile("rgt_sirs_frames_")
  on.exit(unlink(frames_dir, recursive = TRUE), add = TRUE)
  set.seed(19)
  anim_fr <- animate_sirs(g, pos = pos, n_steps = 3L, x = 0.1, r = 0.1, s = 0.05,
                          init_state = init_one, mode = "frames",
                          output = frames_dir, width = 240L, height = 200L)
  expect_equal(length(anim_fr$files), 4L)
  expect_true(all(file.exists(anim_fr$files)))
  # Verify PNG magic bytes on first frame
  sig <- readBin(anim_fr$files[1L], "raw", n = 4L)
  expect_equal(sig, as.raw(c(0x89, 0x50, 0x4e, 0x47)))

  # Invalid input checks
  expect_error(animate_sirs(g, x = 1.5, mode = "none"), "probability")
  expect_error(animate_sirs(g, init_state = rep(3L, n), mode = "none"), "init_state")
})

test_that("animate_dynamic_layout shortens edges under spatial segregation and exports HTML/PNG", {
  # Build a random connected graph with min degree >= 2
  set.seed(314)
  n <- 36L
  g <- graph_create(n, directed = FALSE)
  e_ring <- cbind(seq_len(n), c(seq_len(n - 1L) + 1L, 1L))
  # Add random long-range chords
  chords <- matrix(0L, nrow = 40L, ncol = 2L)
  for (i in seq_len(40L)) {
    u <- sample.int(n, 1L)
    v <- sample.int(n, 1L)
    while (v == u) v <- sample.int(n, 1L)
    chords[i, ] <- c(u, v)
  }
  graph_add_edges(g, rbind(e_ring, chords))

  pos0 <- sfdp_layout(g, K = 0.5, max_iter = 40L)
  html_file <- tempfile("rgt_dyn_", fileext = ".html")
  on.exit(unlink(html_file), add = TRUE)

  res <- animate_dynamic_layout(
    g, pos = pos0, n_steps = 35L, K = 0.5, step = 0.005, n_rewire = 20L,
    mode = "html", output = html_file
  )
  expect_s3_class(res, "rgraphtool_dynamic_layout")
  expect_equal(length(res$positions), 36L)
  expect_equal(length(res$edges), 36L)
  expect_gt(sum(res$rewired), 0L)
  # Distance-shortening rewiring reduces mean Euclidean edge length
  expect_lt(res$mean_edge_length[36L], res$mean_edge_length[1L])

  # Edge count preserved, no self-loops or duplicates, no isolated vertices
  e0_n <- nrow(res$edges[[1L]])
  ef <- res$edges[[36L]]
  expect_equal(nrow(ef), e0_n)
  expect_true(all(ef[, 1L] != ef[, 2L]))
  keys <- paste(pmin(ef[, 1L], ef[, 2L]), pmax(ef[, 1L], ef[, 2L]), sep = "-")
  expect_equal(length(unique(keys)), nrow(ef))
  deg_f <- tabulate(c(ef[, 1L], ef[, 2L]), nbins = n)
  expect_true(all(deg_f >= 1L))

  expect_true(file.exists(html_file))
  txt <- paste(readLines(html_file, warn = FALSE), collapse = "\n")
  expect_match(txt, "meanLen")
  expect_output(print(res), "rgraphtool dynamic layout animation")

  # Frames export
  frames_dir <- tempfile("rgt_dyn_frames_")
  on.exit(unlink(frames_dir, recursive = TRUE), add = TRUE)
  res_fr <- animate_dynamic_layout(
    g, pos = pos0, n_steps = 2L, n_rewire = 5L,
    mode = "frames", output = frames_dir, width = 240L, height = 200L
  )
  expect_equal(length(res_fr$files), 3L)
  expect_true(all(file.exists(res_fr$files)))
})

test_that("interactive_bfs generates interactive HTML and multi-root PNG snapshots", {
  g <- graph_create(15L, directed = FALSE)
  graph_add_edges(g, rbind(
    cbind(1:14, 2:15),
    c(1L, 8L), c(4L, 12L), c(3L, 10L)
  ))
  set.seed(101)
  pos <- sfdp_layout(g, max_iter = 30L)

  html_file <- tempfile("rgt_bfs_", fileext = ".html")
  on.exit(unlink(html_file), add = TRUE)
  ib_html <- interactive_bfs(g, pos = pos, source = 1L, max_dist = 3L,
                             mode = "html", output = html_file)
  expect_s3_class(ib_html, "rgraphtool_interactive_bfs")
  expect_true(file.exists(html_file))
  txt <- paste(readLines(html_file, warn = FALSE), collapse = "\n")
  expect_match(txt, "runBfs")
  expect_match(txt, "mousemove")
  expect_match(txt, "#CE5C00")
  expect_output(print(ib_html), "rgraphtool interactive BFS visualization")

  # Multi-root snapshot to PNG
  snap_file <- tempfile("rgt_bfs_snap_", fileext = ".png")
  on.exit(unlink(snap_file), add = TRUE)
  ib_snap <- interactive_bfs(g, pos = pos, source = c(1L, 5L, 10L), max_dist = 2L,
                             mode = "snapshot", output = snap_file,
                             width = 480L, height = 320L)
  expect_true(file.exists(snap_file))
  expect_gt(file.info(snap_file)$size, 1000)
  expect_equal(ib_snap$bfs$source, 10L)

  # Non-interactive locator mode errors cleanly
  if (!interactive()) {
    expect_error(interactive_bfs(g, pos = pos, mode = "locator"), "interactive R session")
  }
})
