plot_bipartite <- function(state, type = c("matrix", "meso", "network"),
                           level = 1L, weight = NULL,
                           main = NULL,
                           palette = NULL,
                           block_border = "#222222",
                           ...) {
  type <- match.arg(type)
  if (!inherits(state, "rgraphtool_bipartite_block_state") &&
      !inherits(state, "rgraphtool_nested_bipartite_block_state")) {
    stop("state must be a bipartite or nested bipartite block state", call. = FALSE)
  }

  mem <- bipartite_membership(state, level = level)
  b_row <- mem[[1L]]
  b_col <- mem[[2L]]
  dl <- state$dim_labels

  if (type == "matrix") {
    mat <- as_incidence_matrix(state$graph, type = state$type, weight = weight)
    ord_r <- order(b_row)
    ord_c <- order(b_col)
    mat_ord <- mat[ord_r, ord_c, drop = FALSE]
    br_ord <- b_row[ord_r]
    bc_ord <- b_col[ord_c]

    nr <- nrow(mat_ord)
    nc <- ncol(mat_ord)
    if (is.null(palette)) {
      palette <- c("#ffffff", grDevices::hcl.colors(32L, "YlGnBu", rev = TRUE))
    }
    if (is.null(main)) {
      main <- sprintf("Bipartite incidence matrix (%s: %d blocks, %s: %d blocks)",
                      dl[1L], max(b_row, 0L), dl[2L], max(b_col, 0L))
    }

    # Flip rows so row 1 is at the top (standard matrix layout)
    z <- t(mat_ord[rev(seq_len(nr)), , drop = FALSE])
    graphics::image(
      x = seq_len(nc), y = seq_len(nr), z = z,
      col = palette, axes = FALSE,
      xlab = dl[2L], ylab = dl[1L], main = main, ...
    )
    graphics::box()

    # Draw block boundary lines
    if (length(bc_ord) > 1L) {
      cuts_c <- which(diff(bc_ord) != 0L)
      if (length(cuts_c)) graphics::abline(v = cuts_c + 0.5, col = block_border, lwd = 1.4)
    }
    if (length(br_ord) > 1L) {
      cuts_r <- which(diff(rev(br_ord)) != 0L)
      if (length(cuts_r)) graphics::abline(h = cuts_r + 0.5, col = block_border, lwd = 1.4)
    }
    return(invisible(list(matrix = mat_ord, row_order = ord_r, col_order = ord_c)))
  }

  if (type == "meso") {
    dens <- bipartite_block_matrix(state, level = level, normalize = "density")
    counts <- bipartite_edge_counts(state, level = level)
    B0 <- nrow(dens)
    B1 <- ncol(dens)
    n_row <- as.numeric(tabulate(b_row, nbins = B0))
    n_col <- as.numeric(tabulate(b_col, nbins = B1))

    y0 <- if (B0 == 1L) 0.5 else seq(0.9, 0.1, length.out = B0)
    y1 <- if (B1 == 1L) 0.5 else seq(0.9, 0.1, length.out = B1)

    if (is.null(main)) {
      main <- sprintf("Bipartite mesoscale blocks (%s vs %s)", dl[1L], dl[2L])
    }
    graphics::plot.new()
    graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
    graphics::title(main = main, ...)

    max_d <- max(dens, 1e-12)
    for (r in seq_len(B0)) {
      for (s in seq_len(B1)) {
        if (counts[r, s] <= 0) next
        frac <- dens[r, s] / max_d
        lwd <- 0.8 + 5.5 * frac
        col_e <- grDevices::adjustcolor("#3c5488", alpha.f = 0.25 + 0.7 * frac)
        graphics::segments(0.22, y0[r], 0.78, y1[s], lwd = lwd, col = col_e)
      }
    }

    cex0 <- 1.6 + 2.5 * sqrt(n_row / max(n_row, 1))
    cex1 <- 1.6 + 2.5 * sqrt(n_col / max(n_col, 1))
    graphics::points(rep(0.22, B0), y0, pch = 21, bg = "#e64b35", col = "#222222", cex = cex0)
    graphics::points(rep(0.78, B1), y1, pch = 22, bg = "#00a087", col = "#222222", cex = cex1)

    graphics::text(0.14, y0, labels = sprintf("%s %d (n=%d)", dl[1L], seq_len(B0), n_row),
                   adj = c(1, 0.5), cex = 0.85)
    graphics::text(0.86, y1, labels = sprintf("%s %d (n=%d)", dl[2L], seq_len(B1), n_col),
                   adj = c(0, 0.5), cex = 0.85)
    return(invisible(list(density = dens, counts = counts)))
  }

  # type == "network": two-column bipartite node-link layout grouped by block
  rows <- which(state$type == 0L)
  cols <- which(state$type == 1L)
  ord_r <- order(b_row)
  ord_c <- order(b_col)
  pos <- matrix(0, nrow = length(state$type), ncol = 2L)
  pos[rows[ord_r], 1L] <- 0.15
  pos[rows[ord_r], 2L] <- if (length(rows) == 1L) 0.5 else seq(0.95, 0.05, length.out = length(rows))
  pos[cols[ord_c], 1L] <- 0.85
  pos[cols[ord_c], 2L] <- if (length(cols) == 1L) 0.5 else seq(0.95, 0.05, length.out = length(cols))

  B0 <- max(b_row, 1L)
  B1 <- max(b_col, 1L)
  pal0 <- grDevices::hcl.colors(max(B0, 2L), "Warm")[seq_len(B0)]
  pal1 <- grDevices::hcl.colors(max(B1, 2L), "Cold")[seq_len(B1)]
  vcol <- character(length(state$type))
  vcol[rows] <- pal0[b_row]
  vcol[cols] <- pal1[b_col]

  if (is.null(main)) {
    main <- sprintf("Bipartite network (%s vs %s)", dl[1L], dl[2L])
  }
  plot_graph(state$graph, pos = pos, vertex_color = vcol,
             edge_color = grDevices::adjustcolor("#444444", alpha.f = 0.35),
             main = main, ...)
  invisible(pos)
}
