spring_layout <- function(graph, groups = NULL, iter = 120L) {
  n <- graph_num_vertices(graph)
  if (n == 0L) {
    return(matrix(numeric(0), ncol = 2L, dimnames = list(NULL, c("x", "y"))))
  }
  if (n == 1L) {
    return(matrix(c(0, 0), ncol = 2L, dimnames = list(NULL, c("x", "y"))))
  }
  ord <- if (!is.null(groups) && length(groups) == n) {
    order(groups, seq_len(n))
  } else {
    seq_len(n)
  }
  rank <- integer(n)
  rank[ord] <- seq_len(n)
  theta <- 2 * pi * (rank - 1) / n
  pos <- cbind(x = cos(theta), y = sin(theta))
  edges <- graph_edges(graph)
  edges <- edges[edges[, 1] != edges[, 2], , drop = FALSE]
  if (nrow(edges) == 0L || n > 400L) {
    return(pos)
  }
  k <- sqrt(4 / n)
  temp <- 0.3
  cool <- temp / (iter + 1)
  for (step in seq_len(iter)) {
    dx <- outer(pos[, 1], pos[, 1], "-")
    dy <- outer(pos[, 2], pos[, 2], "-")
    dist <- sqrt(dx * dx + dy * dy)
    diag(dist) <- Inf
    dist[dist < 1e-4] <- 1e-4
    rep_f <- (k * k) / (dist * dist)
    disp <- cbind(rowSums(dx * rep_f), rowSums(dy * rep_f))
    u <- edges[, 1]
    v <- edges[, 2]
    ex <- pos[u, 1] - pos[v, 1]
    ey <- pos[u, 2] - pos[v, 2]
    ed <- pmax(sqrt(ex * ex + ey * ey), 1e-4)
    att_f <- ed / k
    ax <- ex * att_f / ed
    ay <- ey * att_f / ed
    for (e in seq_len(nrow(edges))) {
      disp[u[e], 1] <- disp[u[e], 1] - ax[e]
      disp[u[e], 2] <- disp[u[e], 2] - ay[e]
      disp[v[e], 1] <- disp[v[e], 1] + ax[e]
      disp[v[e], 2] <- disp[v[e], 2] + ay[e]
    }
    dlen <- pmax(sqrt(disp[, 1]^2 + disp[, 2]^2), 1e-4)
    step_len <- pmin(dlen, temp)
    pos <- pos + disp * (step_len / dlen)
    temp <- max(temp - cool, 1e-3)
  }
  pos[, 1] <- pos[, 1] - mean(pos[, 1])
  pos[, 2] <- pos[, 2] - mean(pos[, 2])
  sc <- max(abs(pos))
  if (sc > 0) pos <- pos / sc
  colnames(pos) <- c("x", "y")
  pos
}

draw_pie_vertex <- function(x, y, r, probs, cols, border = "grey20", nseg = 36L) {
  nz <- which(probs > 1e-6)
  if (length(nz) == 0L) return(invisible(NULL))
  if (length(nz) == 1L) {
    ang <- seq(0, 2 * pi, length.out = nseg + 1L)
    graphics::polygon(x + r * cos(ang), y + r * sin(ang),
                      col = cols[nz], border = border, lwd = 0.6)
    return(invisible(NULL))
  }
  p <- probs[nz]
  p <- p / sum(p)
  starts <- c(0, cumsum(p)[-length(p)]) * 2 * pi
  ends <- cumsum(p) * 2 * pi
  for (i in seq_along(nz)) {
    k <- max(2L, ceiling(nseg * p[i]))
    ang <- seq(starts[i], ends[i], length.out = k + 1L)
    graphics::polygon(c(x, x + r * cos(ang)), c(y, y + r * sin(ang)),
                      col = cols[nz[i]], border = border, lwd = 0.5)
  }
}

plot_marginals <- function(graph, mode, positions = NULL, radius = 0.04,
                           palette = NULL,
                           edge_col = grDevices::adjustcolor("grey40", 0.35),
                           border = "grey20", labels = NULL, label_cex = 0.65,
                           main = NULL, ...) {
  graph_pointer(graph)
  n <- graph_num_vertices(graph)
  if (inherits(mode, "graphtool_partition_mode")) {
    probs <- mode_marginals(mode, probability = TRUE, level = 0L)
    cons <- if (isTRUE(mode$nested)) mode_consensus(mode)[[1L]] else mode_consensus(mode)
  } else if (is.matrix(mode) && is.numeric(mode)) {
    if (nrow(mode) != n) {
      stop("marginal matrix must have one row per vertex", call. = FALSE)
    }
    rs <- rowSums(mode)
    probs <- mode / pmax(rs, 1)
    probs[rs == 0, ] <- 0
    cons <- max.col(probs, ties.method = "first")
  } else {
    stop("mode must be a partition mode from partition_mode_state() or a numeric matrix",
         call. = FALSE)
  }
  if (nrow(probs) != n) {
    stop("graph and mode must have the same number of vertices", call. = FALSE)
  }
  if (!is.numeric(radius) || length(radius) != 1L || is.na(radius) || radius <= 0) {
    stop("radius must be a positive number", call. = FALSE)
  }
  if (is.null(positions)) {
    positions <- spring_layout(graph, groups = cons)
  } else {
    positions <- as.matrix(positions)
    if (nrow(positions) != n || ncol(positions) != 2L || any(!is.finite(positions))) {
      stop("positions must be a finite numeric matrix with one (x, y) row per vertex",
           call. = FALSE)
    }
    positions[, 1] <- positions[, 1] - mean(positions[, 1])
    positions[, 2] <- positions[, 2] - mean(positions[, 2])
    sc <- max(abs(positions))
    if (sc > 0) positions <- positions / sc
  }
  B <- ncol(probs)
  if (is.null(palette)) {
    palette <- grDevices::hcl.colors(max(B, 1L), "Dynamic")
  } else if (length(palette) < B) {
    palette <- rep_len(palette, B)
  }

  pad <- 1 + 2 * radius
  graphics::plot.new()
  graphics::plot.window(c(-pad, pad), c(-pad, pad), asp = 1, ...)
  edges <- graph_edges(graph)
  edges <- edges[edges[, 1] != edges[, 2], , drop = FALSE]
  if (nrow(edges) > 0L) {
    graphics::segments(positions[edges[, 1], 1], positions[edges[, 1], 2],
                       positions[edges[, 2], 1], positions[edges[, 2], 2],
                       col = edge_col, lwd = 0.8)
  }
  for (v in seq_len(n)) {
    draw_pie_vertex(positions[v, 1], positions[v, 2], radius, probs[v, ],
                    palette, border = border)
  }
  if (!is.null(labels)) {
    graphics::text(positions[, 1], positions[, 2] + radius * 1.3,
                   labels = labels, cex = label_cex)
  }
  if (!is.null(main)) graphics::title(main)
  invisible(list(positions = positions, probabilities = probs, palette = palette))
}
