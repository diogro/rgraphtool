hierarchy_layout <- function(state) {
  validate_nested_state(state)
  levels <- state$levels
  top <- match(1L, vapply(levels, max, integer(1)))
  if (is.na(top)) {
    levels <- c(levels, list(rep(1L, max(levels[[length(levels)]]))))
    top <- length(levels)
  }
  levels <- levels[seq_len(top)]
  n <- length(levels[[1]])
  # ancestors[v, j]: group of vertex v at level j - 1; the last column is the root
  ancestors <- matrix(0L, n, top)
  ancestors[, 1] <- levels[[1]]
  for (j in seq_len(top)[-1]) {
    ancestors[, j] <- levels[[j]][ancestors[, j - 1L]]
  }
  ordering <- do.call(order, c(rev(as.data.frame(ancestors)), list(seq_len(n))))
  rank <- integer(n)
  rank[ordering] <- seq_len(n)
  theta <- 2 * pi * (rank - 1) / n
  vertices <- cbind(x = cos(theta), y = sin(theta))
  groups <- vector("list", top)
  for (j in seq_len(top)) {
    angle <- 2 * pi * (tapply(rank, ancestors[, j], mean) - 1) / n
    radius <- 1 - j / top
    groups[[j]] <- cbind(x = radius * cos(angle), y = radius * sin(angle))
  }
  list(vertices = vertices, groups = groups, ancestors = ancestors, levels = levels)
}

bundle_path <- function(layout, u, v, beta) {
  au <- layout$ancestors[u, ]
  av <- layout$ancestors[v, ]
  lca <- which(au == av)[1]
  up <- lapply(seq_len(lca), function(j) layout$groups[[j]][au[j], ])
  down <- lapply(rev(seq_len(lca - 1L)), function(j) layout$groups[[j]][av[j], ])
  points <- rbind(layout$vertices[u, ], do.call(rbind, up), do.call(rbind, down),
                  layout$vertices[v, ])
  k <- nrow(points)
  t <- (seq_len(k) - 1) / (k - 1)
  straight <- cbind(points[1, 1] + t * (points[k, 1] - points[1, 1]),
                    points[1, 2] + t * (points[k, 2] - points[1, 2]))
  beta * points + (1 - beta) * straight
}

plot_hierarchy <- function(state, beta = 0.8, vertex_col = NULL,
                           edge_col = grDevices::adjustcolor("grey30", 0.25),
                           vertex_cex = 0.6, show_tree = TRUE, main = NULL, ...) {
  layout <- hierarchy_layout(state)
  if (!is.numeric(beta) || length(beta) != 1L || is.na(beta) || beta < 0 || beta > 1) {
    stop("beta must be a number between 0 and 1", call. = FALSE)
  }
  groups0 <- layout$levels[[1]]
  if (is.null(vertex_col)) {
    vertex_col <- grDevices::hcl.colors(max(groups0), "Dynamic")[groups0]
  }
  graphics::plot.new()
  graphics::plot.window(c(-1.05, 1.05), c(-1.05, 1.05), asp = 1)
  edges <- graph_edges(state$graph)
  edges <- edges[edges[, 1] != edges[, 2], , drop = FALSE]
  edge_col <- rep_len(edge_col, nrow(edges))
  for (i in seq_len(nrow(edges))) {
    path <- bundle_path(layout, edges[i, 1], edges[i, 2], beta)
    graphics::xspline(path[, 1], path[, 2], shape = 1, open = TRUE, border = edge_col[i])
  }
  if (show_tree) {
    depth <- length(layout$groups)
    for (j in seq_len(depth)) {
      child <- if (j == 1L) layout$vertices else layout$groups[[j - 1L]]
      parent_index <- if (j == 1L) groups0 else layout$levels[[j]]
      parent <- layout$groups[[j]][parent_index, , drop = FALSE]
      graphics::segments(child[, 1], child[, 2], parent[, 1], parent[, 2],
                         col = grDevices::adjustcolor("firebrick", 0.5), lwd = 0.8)
      graphics::points(layout$groups[[j]], pch = 15, cex = 0.7, col = "firebrick")
    }
  }
  graphics::points(layout$vertices, pch = 19, cex = vertex_cex, col = vertex_col, ...)
  if (!is.null(main)) graphics::title(main)
  invisible(layout)
}
