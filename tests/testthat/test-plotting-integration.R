test_that("prop_to_size and map_colors scale properties accurately", {
  x <- c(1, 4, 9, 16)
  s <- prop_to_size(x, mi = 1, ma = 4, power = 0.5)
  expect_equal(min(s), 1)
  expect_equal(max(s), 4)
  expect_true(all(diff(s) > 0))

  # Log scaling
  s_log <- prop_to_size(c(1, 10, 100), mi = 1, ma = 3, power = 1, log = TRUE)
  expect_equal(s_log, c(1, 2, 3), tolerance = 1e-12)
  expect_error(prop_to_size(c(0, 1, 2), log = TRUE), "strictly positive")

  # Constant and empty vectors
  expect_equal(prop_to_size(c(5, 5, 5), mi = 1, ma = 3), c(2, 2, 2))
  expect_equal(prop_to_size(numeric(0)), numeric(0))

  # Color mapping + limits clipping
  cols <- map_colors(c(0, 5, 10, 20), palette = "Viridis", limits = c(0, 10))
  expect_s3_class(cols, "rgraphtool_color_map")
  expect_equal(length(cols), 4L)
  expect_equal(cols[3L], cols[4L])
  expect_equal(attr(cols, "limits"), c(0, 10))
})

test_that("save_graph_plot exports vector (PDF, SVG, EPS) and raster (PNG, JPEG, TIFF) formats", {
  g <- graph_create(10L, directed = TRUE)
  graph_add_edges(g, rbind(cbind(1:9, 2:10), c(10L, 1L), c(1L, 5L)))
  set.seed(11)
  pos <- sfdp_layout(g, max_iter = 30L)
  deg <- graph_degree(g)
  vcol <- map_colors(deg, palette = "Plasma")
  vsiz <- prop_to_size(deg, mi = 0.9, ma = 2.0)

  # PDF vector output
  f_pdf <- tempfile("rgt_vec_", fileext = ".pdf")
  on.exit(unlink(f_pdf), add = TRUE)
  plot_graph(g, pos = pos, vertex_color = vcol, vertex_size = vsiz,
             edge_curve = 0.12, axes = TRUE, grid_lines = TRUE,
             xlab = "x", ylab = "y", vorder = deg,
             main = "Vector PDF", output = f_pdf, width = 6, height = 4)
  expect_true(file.exists(f_pdf))
  expect_equal(rawToChar(readBin(f_pdf, "raw", n = 4L)), "%PDF")

  # EPS vector output
  f_eps <- tempfile("rgt_vec_", fileext = ".eps")
  on.exit(unlink(f_eps), add = TRUE)
  save_graph_plot(plot_colorbar(vcol, title = "Degree"), output = f_eps, width = 2, height = 4)
  expect_true(file.exists(f_eps))
  expect_equal(rawToChar(readBin(f_eps, "raw", n = 4L)), "%!PS")

  # SVG vector output (via svglite or grDevices::svg)
  if (requireNamespace("svglite", quietly = TRUE)) {
    f_svg <- tempfile("rgt_vec_", fileext = ".svg")
    on.exit(unlink(f_svg), add = TRUE)
    plot_graph(g, pos = pos, vertex_color = vcol, edge_curve = 0.1,
               output = f_svg, width = 5, height = 4)
    expect_true(file.exists(f_svg))
    svg_head <- paste(readLines(f_svg, n = 10L, warn = FALSE), collapse = "\n")
    expect_match(svg_head, "<svg")
  }

  # PNG raster output
  f_png <- tempfile("rgt_ras_", fileext = ".png")
  on.exit(unlink(f_png), add = TRUE)
  save_graph_plot({
    op <- graphics::par(mfrow = c(1L, 2L))
    on.exit(graphics::par(op), add = TRUE)
    plot_graph(g, pos = pos, vertex_color = vcol, axes = TRUE, main = "Panel A")
    plot_colorbar(vcol, horizontal = TRUE, title = "Degree")
  }, output = f_png, width = 7, height = 3.5, dpi = 100)
  expect_true(file.exists(f_png))
  expect_equal(readBin(f_png, "raw", n = 4L), as.raw(c(0x89, 0x50, 0x4e, 0x47)))

  # JPEG and TIFF raster outputs
  f_jpg <- tempfile("rgt_ras_", fileext = ".jpg")
  f_tif <- tempfile("rgt_ras_", fileext = ".tiff")
  on.exit(unlink(c(f_jpg, f_tif)), add = TRUE)
  plot_graph(g, pos = pos, output = f_jpg, width = 4, height = 3, dpi = 72)
  plot_graph(g, pos = pos, output = f_tif, width = 4, height = 3, dpi = 72)
  expect_equal(readBin(f_jpg, "raw", n = 2L), as.raw(c(0xff, 0xd8)))
  expect_gt(file.info(f_tif)$size, 500)

  # Invalid extension
  expect_error(save_graph_plot(NULL, output = tempfile(fileext = ".xyz")), "Unsupported")
})

test_that("graph_grob, graph_to_data_frames, and ggplot_graph integrate with grid and ggplot2", {
  g <- graph_create(8L, directed = TRUE)
  graph_add_edges(g, rbind(cbind(1:7, 2:8), c(8L, 1L), c(2L, 6L)),
                  properties = list(weight = seq_len(9L)))
  set.seed(5)
  pos <- sfdp_layout(g, max_iter = 25L)

  # Tidy data frames
  dfs <- graph_to_data_frames(g, pos = pos, vertex_data = list(group = rep(c("A", "B"), each = 4L)))
  expect_equal(nrow(dfs$vertices), 8L)
  expect_equal(nrow(dfs$edges), 9L)
  expect_true(all(c("id", "x", "y", "degree", "group") %in% names(dfs$vertices)))
  expect_true(all(c("from", "to", "x", "y", "xend", "yend", "weight") %in% names(dfs$edges)))

  # grid grob rendering
  gb <- graph_grob(g, pos = pos, vertex_color = map_colors(seq_len(8L)),
                   vertex_halo = c(TRUE, rep(FALSE, 7L)),
                   edge_curve = 0.15, axes = TRUE, bg = "#F8FAFC", main = "Grid Grob")
  expect_s3_class(gb, "rgraphtool_grob")
  expect_s3_class(gb, "grob")

  f_grid <- tempfile("rgt_grid_", fileext = ".pdf")
  on.exit(unlink(f_grid), add = TRUE)
  save_graph_plot(gb, output = f_grid, width = 5, height = 4)
  expect_true(file.exists(f_grid))
  expect_gt(file.info(f_grid)$size, 500)

  # ggplot2 integration (if installed)
  skip_if_not_installed("ggplot2")
  p1 <- ggplot_graph(g, pos = pos, vertex_color = graph_degree(g),
                     vertex_size = graph_degree(g), edge_curve = 0.15,
                     axes = TRUE, main = "ggplot2 Numeric", color_label = "Degree")
  expect_s3_class(p1, "ggplot")

  p2 <- ggplot_graph(g, pos = pos, vertex_color = factor(rep(c("C1", "C2"), each = 4L)),
                     vertex_label = paste0("V", seq_len(8L)),
                     edge_color = seq_len(9L), main = "ggplot2 Factor")
  expect_s3_class(p2, "ggplot")

  f_gg <- tempfile("rgt_gg_", fileext = ".png")
  on.exit(unlink(f_gg), add = TRUE)
  save_graph_plot(p1, output = f_gg, width = 5, height = 4, dpi = 96)
  expect_true(file.exists(f_gg))
  expect_gt(file.info(f_gg)$size, 1000)
})

test_that("project_coords, great_circle_arc, world_coastlines, and plot_graph_map work across projections", {
  # Projection math checks
  lons <- c(-120, 0, 120, 180)
  lats <- c(-30, 0, 45, 0)
  for (proj in c("robinson", "mercator", "equirectangular", "mollweide", "orthographic")) {
    pr <- project_coords(lons, lats, projection = proj, lon0 = 0, lat0 = 0)
    expect_equal(nrow(pr), 4L)
    expect_true(all(c("x", "y", "visible") %in% names(pr)))
    expect_true(all(is.finite(pr$x)))
    expect_true(all(is.finite(pr$y)))
  }

  # Orthographic hemisphere visibility check
  pr_ortho <- project_coords(c(0, 180), c(0, 0), projection = "orthographic", lon0 = 0, lat0 = 0)
  expect_equal(pr_ortho$visible, c(TRUE, FALSE))

  # Mollweide exact bounds
  pr_moll <- project_coords(c(-180, 0, 180, 0, 0), c(0, 0, 0, 90, -90), projection = "mollweide")
  expect_equal(pr_moll$y[c(4L, 5L)], c(sqrt(2), -sqrt(2)), tolerance = 1e-6)
  expect_equal(abs(pr_moll$x[1L]), 2 * sqrt(2), tolerance = 1e-6)

  # Great-circle arc across Pacific splits across antimeridian cleanly
  arc <- great_circle_arc(139.7, 35.7, -118.2, 34.0, n = 30L)
  expect_equal(nrow(arc), 30L)
  expect_equal(arc$lon[1L], 139.7, tolerance = 1e-6)
  expect_equal(arc$lat[1L], 35.7, tolerance = 1e-6)
  expect_equal(arc$lon[30L], -118.2, tolerance = 1e-6)
  expect_equal(arc$lat[30L], 34.0, tolerance = 1e-6)
  expect_gt(length(unique(arc$segment)), 1L)
  for (sid in unique(arc$segment)) {
    sub <- arc[arc$segment == sid, , drop = FALSE]
    if (nrow(sub) >= 2L) expect_true(all(abs(diff(sub$lon)) <= 180))
  }

  # Built-in offline coastlines
  wc <- world_coastlines()
  expect_true(is.data.frame(wc))
  expect_true(all(c("lon", "lat", "group", "region") %in% names(wc)))
  expect_gt(length(unique(wc$region)), 6L)

  # Full global & regional map rendering (PDF + PNG)
  cities_lon <- c(-74.0, -0.1, 139.7, -118.2, -46.6, 151.2, 2.35, 28.0)
  cities_lat <- c(40.7, 51.5, 35.7, 34.0, -23.5, -33.9, 48.9, -26.2)
  cities_lbl <- c("NYC", "LON", "TYO", "LAX", "GRU", "SYD", "PAR", "JNB")
  g_air <- graph_create(length(cities_lon), directed = FALSE)
  graph_add_edges(g_air, rbind(
    c(1L, 2L), c(1L, 4L), c(1L, 5L), c(2L, 7L), c(2L, 8L),
    c(2L, 3L), c(3L, 4L), c(3L, 6L), c(5L, 8L), c(4L, 6L)
  ))

  f_map_pdf <- tempfile("rgt_map_", fileext = ".pdf")
  f_map_png <- tempfile("rgt_map_", fileext = ".png")
  on.exit(unlink(c(f_map_pdf, f_map_png)), add = TRUE)

  res_rob <- plot_graph_map(
    g_air, lon = cities_lon, lat = cities_lat, projection = "robinson",
    vertex_size = prop_to_size(graph_degree(g_air), 1.0, 2.2),
    vertex_label = cities_lbl, main = "Global Hub Network (Robinson)",
    output = f_map_pdf, width = 8, height = 4.5
  )
  expect_true(file.exists(f_map_pdf))
  expect_equal(rawToChar(readBin(f_map_pdf, "raw", n = 4L)), "%PDF")
  expect_equal(nrow(res_rob$projected), length(cities_lon))

  res_ortho <- plot_graph_map(
    g_air, lon = cities_lon, lat = cities_lat, projection = "orthographic",
    lon0 = -20, lat0 = 25, geodesic = TRUE,
    main = "Orthographic Globe", output = f_map_png, width = 5, height = 5, dpi = 96
  )
  expect_true(file.exists(f_map_png))
  expect_equal(readBin(f_map_png, "raw", n = 4L), as.raw(c(0x89, 0x50, 0x4e, 0x47)))

  # Regional bounding box (North Atlantic / Europe)
  f_reg <- tempfile("rgt_map_reg_", fileext = ".png")
  on.exit(unlink(f_reg), add = TRUE)
  plot_graph_map(
    g_air, lon = cities_lon, lat = cities_lat, projection = "mercator",
    xlim = c(-90, 30), ylim = c(20, 65), geodesic = FALSE,
    output = f_reg, width = 6, height = 4, dpi = 72
  )
  expect_true(file.exists(f_reg))

  # Validation error checks
  expect_error(plot_graph_map(g_air, lon = 1:3, lat = 1:3), "length")
})

test_that("plot_hierarchy marks groups with arcs outside the circle", {
  g <- graph_create(12, directed = FALSE)
  graph_add_edges(g, do.call(rbind, lapply(list(1:4, 5:8, 9:12), function(v) t(combn(v, 2)))))
  state <- nested_block_state(g, list(rep(1:3, each = 4), c(1, 1, 2), c(1, 1)))
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  a0 <- plot_hierarchy(state, arc_level = 0)$arcs
  expect_equal(a0$size, c(4L, 4L, 4L))
  expect_equal(a0$label, c("1", "2", "3"))
  # arcs do not overlap and stay inside one turn
  o <- order(a0$start)
  expect_true(all(a0$end[o][-3] < a0$start[o][-1]))
  expect_true(all(a0$start >= -pi / 12 & a0$end <= 2 * pi))
  a1 <- plot_hierarchy(state, arc_level = 1, arc_labels = letters[1:12])$arcs
  expect_equal(a1$size, c(8L, 4L))
  expect_equal(a1$label, c("a", "i"))
  expect_equal(plot_hierarchy(state, arc_level = 1, arc_col = c("1" = "red", "2" = "blue"))$arcs$col,
               c("red", "blue"))
  expect_error(plot_hierarchy(state, arc_level = 3), "arc_level")
  expect_error(plot_hierarchy(state, arc_level = 0, arc_labels = 1:3), "one label per vertex")
})
