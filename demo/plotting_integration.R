# Cookbook demo: Multi-panel plotting, grid/ggplot2 integration, and vector/raster export
# R-native counterpart of graph-tool's matplotlib integration examples.
# Runs fully offline without requiring Python.

library(rgraphtool)

set.seed(2026)

# Construct a 3-community modular network (60 vertices) and fit an SBM
n_blocks <- 3L
b_size <- 20L
n <- n_blocks * b_size
true_b <- rep(seq_len(n_blocks), each = b_size)

edge_list <- list()
for (i in seq_len(n - 1L)) {
  for (j in (i + 1L):n) {
    p_edge <- if (true_b[i] == true_b[j]) 0.25 else 0.02
    if (stats::runif(1L) < p_edge) {
      edge_list[[length(edge_list) + 1L]] <- c(i, j)
    }
  }
}
edges <- do.call(rbind, edge_list)
g <- graph_create(n, directed = FALSE)
graph_add_edges(g, edges)

fit <- minimize_blockmodel_dl(g, degree_corrected = TRUE)
groups <- block_membership(fit)
pos <- sfdp_layout(g, groups = groups, mu = 0.4, max_iter = 100L)

deg <- graph_degree(g)
vcol_deg <- map_colors(deg, palette = "Plasma")
vsiz_deg <- prop_to_size(deg, mi = 0.85, ma = 2.1, power = 0.5)
vcol_grp <- grDevices::hcl.colors(max(groups), palette = "Dark 3")[groups]

# 1. Multi-panel base-R figure (graph + colorbar + SBM communities + degree histogram)
#    exported to both vector PDF and raster PNG.
draw_multipanel <- function() {
  graphics::layout(matrix(c(1, 2, 3, 4), nrow = 2L, ncol = 2L, byrow = TRUE),
                   widths = c(4, 1.3))
  plot_graph(
    g, pos = pos,
    vertex_color = vcol_deg,
    vertex_size = vsiz_deg,
    vorder = deg,
    edge_curve = 0.12,
    axes = TRUE,
    grid_lines = TRUE,
    xlab = "Layout x",
    ylab = "Layout y",
    main = "Degree-Scaled Layout with Axes & Grid"
  )
  plot_colorbar(vcol_deg, title = "Degree")
  plot_graph(
    g, pos = pos,
    vertex_color = vcol_grp,
    vertex_size = vsiz_deg,
    edge_curve = 0.08,
    axes = TRUE,
    xlab = "Layout x",
    ylab = "Layout y",
    main = sprintf("Inferred SBM (B = %d, DL = %.1f)", max(groups),
                   block_description_length(fit))
  )
  op <- graphics::par(mar = c(4.0, 4.0, 2.5, 1.2))
  on.exit(graphics::par(op), add = TRUE)
  graphics::hist(deg, breaks = 10L, col = "#3B82F6", border = "white",
                 main = "Degree Distribution", xlab = "Degree", ylab = "Count")
}

out_pdf <- file.path(tempdir(), "rgraphtool_plotting_multipanel.pdf")
out_png <- file.path(tempdir(), "rgraphtool_plotting_multipanel.png")
save_graph_plot(draw_multipanel(), output = out_pdf, width = 9, height = 7)
save_graph_plot(draw_multipanel(), output = out_png, width = 9, height = 7, dpi = 150)

cat(sprintf("Wrote vector PDF (%d bytes): %s\n", file.info(out_pdf)$size, out_pdf))
cat(sprintf("Wrote raster PNG (%d bytes): %s\n", file.info(out_png)$size, out_png))

# 2. Vector SVG output (when svglite is available)
if (requireNamespace("svglite", quietly = TRUE)) {
  out_svg <- file.path(tempdir(), "rgraphtool_plotting_graph.svg")
  plot_graph(g, pos = pos, vertex_color = vcol_deg, vertex_size = vsiz_deg,
             edge_curve = 0.12, axes = TRUE, main = "SVG Vector Export",
             output = out_svg, width = 6, height = 5)
  cat(sprintf("Wrote vector SVG (%d bytes): %s\n", file.info(out_svg)$size, out_svg))
}

# 3. grid grob embedding (two graph_grob panels composed with grid viewports)
g1 <- graph_grob(g, pos = pos, vertex_color = vcol_deg, vertex_size = vsiz_deg,
                 edge_curve = 0.1, axes = TRUE, bg = "#F8FAFC", main = "Grid: Degree")
g2 <- graph_grob(g, pos = pos, vertex_color = vcol_grp, vertex_size = vsiz_deg,
                 axes = TRUE, bg = "#F8FAFC", main = "Grid: SBM Blocks")
out_grid_pdf <- file.path(tempdir(), "rgraphtool_plotting_grid.pdf")
save_graph_plot({
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1L, 2L)))
  grid::pushViewport(grid::viewport(layout.pos.row = 1L, layout.pos.col = 1L))
  grid::grid.draw(g1)
  grid::popViewport()
  grid::pushViewport(grid::viewport(layout.pos.row = 1L, layout.pos.col = 2L))
  grid::grid.draw(g2)
  grid::popViewport(2L)
}, output = out_grid_pdf, width = 10, height = 4.8)
cat(sprintf("Wrote grid-composed PDF (%d bytes): %s\n",
            file.info(out_grid_pdf)$size, out_grid_pdf))

# 4. Optional ggplot2 integration
if (requireNamespace("ggplot2", quietly = TRUE)) {
  p_gg <- ggplot_graph(
    g, pos = pos,
    vertex_color = deg,
    vertex_size = deg,
    edge_curve = 0.12,
    axes = TRUE,
    palette = "Viridis",
    color_label = "Degree",
    size_label = "Degree",
    main = "rgraphtool + ggplot2 Integration",
    subtitle = sprintf("60-node modular network (%d edges)", graph_num_edges(g))
  )
  out_gg <- file.path(tempdir(), "rgraphtool_plotting_ggplot2.png")
  save_graph_plot(p_gg, output = out_gg, width = 7, height = 5.5, dpi = 150)
  cat(sprintf("Wrote ggplot2 PNG (%d bytes): %s\n", file.info(out_gg)$size, out_gg))
}
