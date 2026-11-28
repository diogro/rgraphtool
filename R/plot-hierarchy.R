hierarchy_layout <- function(state) {
  if (inherits(state, "rgraphtool_nested_weighted_block_state")) {
    validate_nested_weighted_state(state)
  } else if (inherits(state, "rgraphtool_nested_layered_block_state")) {
    validate_nested_layered_state(state)
  } else if (inherits(state, "rgraphtool_nested_ranked_block_state")) {
    validate_nested_ranked_state(state)
  } else {
    validate_nested_state(state)
  }
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

# Arcs outside the vertex circle, one per group at hierarchy level `level`
# (0 = the groups of the vertices). The hierarchy ordering makes every group's
# vertices contiguous, so each group spans one interval of angles. `gap` is the
# fraction of a vertex slot left empty at each end of an arc.
hierarchy_arcs <- function(layout, level, gap = 0.3) {
  depth <- ncol(layout$ancestors)
  if (!is.numeric(level) || length(level) != 1L || is.na(level) || level < 0 ||
      level != floor(level) || level > depth - 1L) {
    stop(sprintf("arc_level must be an integer between 0 and %d", depth - 1L), call. = FALSE)
  }
  n <- nrow(layout$vertices)
  slot <- 2 * pi / n
  rank <- (atan2(layout$vertices[, 2], layout$vertices[, 1]) %% (2 * pi)) / slot
  rank <- round(rank)  # 0-based position on the circle
  grp <- layout$ancestors[, level + 1L]
  ids <- sort(unique(grp))
  lo <- vapply(ids, function(g) min(rank[grp == g]), numeric(1))
  hi <- vapply(ids, function(g) max(rank[grp == g]), numeric(1))
  pad <- min(max(gap, 0), 0.49) * slot
  data.frame(group = ids, size = as.integer(tabulate(match(grp, ids), length(ids))),
             start = (lo - 0.5) * slot + pad, end = (hi + 0.5) * slot - pad,
             mid = (lo + hi) / 2 * slot)
}

#' @noRd
#' @export
plot_hierarchy <- function(state, beta = 0.8, vertex_col = NULL,
                           edge_col = grDevices::adjustcolor("grey30", 0.25),
                           vertex_cex = 0.6, show_tree = TRUE, main = NULL,
                           tree_col = "#729fcf", edge_lwd = 1,
                           arc_level = NULL, arc_col = NULL, arc_labels = TRUE,
                           arc_radius = 1.06, arc_lwd = 5, arc_label_cex = 0.8,
                           arc_label_radius = arc_radius + 0.05, ...) {
  layout <- hierarchy_layout(state)
  if (!is.numeric(beta) || length(beta) != 1L || is.na(beta) || beta < 0 || beta > 1) {
    stop("beta must be a number between 0 and 1", call. = FALSE)
  }
  groups0 <- layout$levels[[1]]
  default_vertex_col <- is.null(vertex_col)
  if (default_vertex_col) {
    vertex_col <- grDevices::hcl.colors(max(groups0), "Dynamic")[groups0]
  }
  arcs <- NULL
  lim <- 1.05
  if (!is.null(arc_level)) {
    arcs <- hierarchy_arcs(layout, arc_level)
    n_arc <- nrow(arcs)
    if (is.null(arc_col)) {
      # level 0 with default vertex colours: the arcs repeat the vertex colours;
      # coarser levels alternate two greys so they are not read as level-0 groups
      arc_col <- if (arc_level == 0 && default_vertex_col) {
        grDevices::hcl.colors(max(groups0), "Dynamic")[arcs$group]
      } else {
        rep_len(c("#4B5563", "#9CA3AF"), n_arc)
      }
    } else if (!is.null(names(arc_col)) && all(as.character(arcs$group) %in% names(arc_col))) {
      arc_col <- arc_col[as.character(arcs$group)]
    } else {
      arc_col <- rep_len(arc_col, n_arc)
    }
    arcs$col <- unname(as.character(arc_col))
    # labels: TRUE = group ids; a vector over vertices labels each arc with the
    # value at its vertices (the first, when they differ)
    lab <- NULL
    if (isTRUE(arc_labels)) {
      lab <- as.character(arcs$group)
    } else if (!isFALSE(arc_labels) && !is.null(arc_labels)) {
      if (length(arc_labels) != nrow(layout$vertices)) {
        stop("arc_labels must be TRUE, FALSE, or one label per vertex", call. = FALSE)
      }
      grp <- layout$ancestors[, arc_level + 1L]
      lab <- as.character(arc_labels)[match(arcs$group, grp)]
    }
    arcs$label <- if (is.null(lab)) NA_character_ else lab
    lim <- max(lim, arc_radius + 0.05, if (!is.null(lab)) arc_label_radius + 0.2)
  }
  graphics::plot.new()
  graphics::plot.window(c(-lim, lim), c(-lim, lim), asp = 1)
  edges <- graph_edges(state$graph)
  edges <- edges[edges[, 1] != edges[, 2], , drop = FALSE]
  # colours and widths follow the non-self-loop edges in graph order, recycled
  edge_col <- rep_len(edge_col, nrow(edges))
  edge_lwd <- rep_len(as.numeric(edge_lwd), nrow(edges))
  # draw thin edges first so the heavy ones stay on top
  for (i in order(edge_lwd)) {
    path <- bundle_path(layout, edges[i, 1], edges[i, 2], beta)
    graphics::xspline(path[, 1], path[, 2], shape = 1, open = TRUE, border = edge_col[i],
                      lwd = edge_lwd[i])
  }
  if (show_tree) {
    depth <- length(layout$groups)
    for (j in seq_len(depth)) {
      child <- if (j == 1L) layout$vertices else layout$groups[[j - 1L]]
      parent_index <- if (j == 1L) groups0 else layout$levels[[j]]
      parent <- layout$groups[[j]][parent_index, , drop = FALSE]
      graphics::segments(child[, 1], child[, 2], parent[, 1], parent[, 2],
                         col = grDevices::adjustcolor(tree_col, 0.6), lwd = 0.8)
      graphics::points(layout$groups[[j]], pch = 15, cex = 0.7,
                       col = grDevices::adjustcolor(tree_col, 0.6))
    }
  }
  graphics::points(layout$vertices, pch = 19, cex = vertex_cex, col = vertex_col, ...)
  if (!is.null(arcs)) {
    for (k in seq_len(nrow(arcs))) {
      a <- seq(arcs$start[k], arcs$end[k], length.out = max(2L, ceiling((arcs$end[k] - arcs$start[k]) * 60)))
      graphics::lines(arc_radius * cos(a), arc_radius * sin(a), col = arcs$col[k],
                      lwd = arc_lwd, lend = "butt", xpd = NA)
      if (!is.na(arcs$label[k])) {
        m <- arcs$mid[k]
        graphics::text(arc_label_radius * cos(m), arc_label_radius * sin(m), arcs$label[k],
                       adj = c((1 - cos(m)) / 2, (1 - sin(m)) / 2), cex = arc_label_cex,
                       xpd = NA, col = "#1F2937")
      }
    }
    layout$arcs <- arcs
  }
  if (!is.null(main)) graphics::title(main)
  invisible(layout)
}

#' @noRd
#' @export
plot_ranked_hierarchy <- function(state,
                                  upstream_col = grDevices::adjustcolor("#4878CF", 0.35),
                                  downstream_col = grDevices::adjustcolor("#D65F5F", 0.55),
                                  lateral_col = grDevices::adjustcolor("#333355", 0.25),
                                  vertex_col = NULL, vertex_cex = 0.85,
                                  curve = 0.15, show_ranks = TRUE,
                                  show_legend = TRUE, main = NULL, ...) {
  b <- ranked_base_blocks(state)
  n <- length(b)
  if (n == 0L) stop("graph has no vertices to plot", call. = FALSE)
  u_v <- state$u[b]
  uniq_u <- sort(unique(u_v))
  k_ranks <- length(uniq_u)
  rank_idx <- match(u_v, uniq_u)
  y_pos <- if (k_ranks > 1L) (rank_idx - 1) / (k_ranks - 1) else rep(0.5, n)

  x_pos <- numeric(n)
  for (rk in seq_len(k_ranks)) {
    idx <- which(rank_idx == rk)
    ord <- idx[order(b[idx], idx)]
    m <- length(ord)
    if (m == 1L) {
      x_pos[ord] <- 0
    } else {
      span <- min(0.92, 0.14 * m)
      x_pos[ord] <- seq(-span, span, length.out = m)
    }
  }
  coords <- cbind(x = x_pos, y = y_pos)

  if (is.null(vertex_col)) {
    pal <- grDevices::hcl.colors(max(b, 1L), "Dynamic")
    vertex_col <- pal[b]
  }

  graphics::plot.new()
  graphics::plot.window(c(-1.05, 1.05), c(-0.08, 1.10))

  if (isTRUE(show_ranks)) {
    y_tiers <- if (k_ranks > 1L) seq(0, 1, length.out = k_ranks) else 0.5
    graphics::abline(h = y_tiers, col = grDevices::adjustcolor("grey75", 0.5), lty = 3)
  }

  ed <- graph_edges(state$graph)
  if (nrow(ed) > 0L) {
    keep <- ed[, 1L] != ed[, 2L]
    ed <- ed[keep, , drop = FALSE]
  }
  if (nrow(ed) > 0L) {
    us <- u_v[ed[, 1L]]
    ut <- u_v[ed[, 2L]]
    dirs <- ifelse(us < ut, 1L, ifelse(us > ut, -1L, 0L))
    draw_order <- c(which(dirs == 0L), which(dirs == 1L), which(dirs == -1L))
    for (i in draw_order) {
      p1 <- coords[ed[i, 1L], ]
      p2 <- coords[ed[i, 2L], ]
      d <- dirs[i]
      col_e <- if (d == 1L) upstream_col else if (d == -1L) downstream_col else lateral_col
      dx <- p2[1L] - p1[1L]
      dy <- p2[2L] - p1[2L]
      dist_e <- sqrt(dx * dx + dy * dy)
      if (dist_e > 1e-8 && curve != 0) {
        sgn <- if (d == -1L) -1 else 1
        mx <- 0.5 * (p1[1L] + p2[1L]) - sgn * curve * dy
        my <- 0.5 * (p1[2L] + p2[2L]) + sgn * curve * dx
        graphics::xspline(c(p1[1L], mx, p2[1L]), c(p1[2L], my, p2[2L]),
                          shape = 1, open = TRUE, border = col_e)
      } else {
        graphics::segments(p1[1L], p1[2L], p2[1L], p2[2L], col = col_e)
      }
    }
  }

  graphics::points(coords, pch = 21, cex = vertex_cex, bg = vertex_col, col = "grey20", ...)

  if (isTRUE(show_legend)) {
    ec <- ranked_edge_counts(state)
    graphics::legend(
      "top", horiz = TRUE, bty = "n", cex = 0.75,
      legend = c(sprintf("upstream (%d)", ec$upstream),
                 sprintf("lateral (%d)", ec$lateral),
                 sprintf("downstream (%d)", ec$downstream)),
      col = c("#4878CF", "#333355", "#D65F5F"), lwd = 2
    )
  }
  if (!is.null(main)) graphics::title(main = main)
  invisible(list(vertices = coords, blocks = b, u = state$u, ranks = u_v))
}
