# Plotting integration: property scaling, colorbars, vector/raster export,
# grid grob embedding, tidy data frames, and ggplot2 integration.

#' @noRd
#' @export
prop_to_size <- function(x, mi = 0.6, ma = 2.4, power = 0.5, log = FALSE) {
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x))) {
    stop("x must be a finite numeric vector", call. = FALSE)
  }
  if (!is.numeric(mi) || length(mi) != 1L || !is.finite(mi) ||
      !is.numeric(ma) || length(ma) != 1L || !is.finite(ma) || mi > ma) {
    stop("mi and ma must be finite scalars with mi <= ma", call. = FALSE)
  }
  power <- validate_pos_scalar(power, "power")
  if (!is.logical(log) || length(log) != 1L || is.na(log)) {
    stop("log must be TRUE or FALSE", call. = FALSE)
  }
  if (length(x) == 0L) return(numeric(0))
  vals <- as.numeric(x)
  if (isTRUE(log)) {
    if (any(vals <= 0)) {
      stop("all values of x must be strictly positive when log = TRUE", call. = FALSE)
    }
    vals <- base::log(vals)
  }
  rng <- range(vals)
  span <- rng[2L] - rng[1L]
  if (span <= 0) {
    return(rep(0.5 * (mi + ma), length(vals)))
  }
  u <- ((vals - rng[1L]) / span)^power
  mi + u * (ma - mi)
}

#' @noRd
#' @export
map_colors <- function(x, palette = "Viridis", n = 256L, limits = NULL,
                       alpha = 1, rev = FALSE) {
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x))) {
    stop("x must be a finite numeric vector", call. = FALSE)
  }
  n <- validate_nonneg_int(n, "n", 2L)
  alpha <- validate_prob_scalar(alpha, "alpha")
  if (!is.logical(rev) || length(rev) != 1L || is.na(rev)) {
    stop("rev must be TRUE or FALSE", call. = FALSE)
  }
  if (length(x) == 0L) {
    out <- character(0)
    attr(out, "limits") <- if (is.null(limits)) c(0, 1) else limits
    attr(out, "palette") <- palette
    attr(out, "alpha") <- alpha
    attr(out, "rev") <- rev
    class(out) <- "rgraphtool_color_map"
    return(out)
  }
  if (is.null(limits)) {
    limits <- range(x)
  } else if (!is.numeric(limits) || length(limits) != 2L ||
             any(!is.finite(limits)) || limits[1L] > limits[2L]) {
    stop("limits must be a finite numeric vector c(min, max) with min <= max", call. = FALSE)
  }
  cols <- grDevices::hcl.colors(n, palette = palette, alpha = alpha, rev = rev)
  span <- limits[2L] - limits[1L]
  if (span <= 0) {
    idx <- rep(as.integer(ceiling(n / 2)), length(x))
  } else {
    u <- (x - limits[1L]) / span
    u <- pmin(1, pmax(0, u))
    idx <- as.integer(floor(u * (n - 1L))) + 1L
  }
  out <- cols[idx]
  attr(out, "limits") <- limits
  attr(out, "palette") <- palette
  attr(out, "alpha") <- alpha
  attr(out, "rev") <- rev
  class(out) <- "rgraphtool_color_map"
  out
}

#' @noRd
#' @export
plot_colorbar <- function(limits, palette = "Viridis", n = 100L, title = NULL,
                          alpha = 1, rev = FALSE, horizontal = FALSE,
                          n_ticks = 5L, digits = 2L, cex_axis = 0.85, ...) {
  if (inherits(limits, "rgraphtool_color_map")) {
    palette <- attr(limits, "palette")
    alpha <- attr(limits, "alpha")
    rev <- attr(limits, "rev")
    limits <- attr(limits, "limits")
  }
  if (!is.numeric(limits) || length(limits) != 2L ||
      any(!is.finite(limits)) || limits[1L] > limits[2L]) {
    stop("limits must be c(min, max) or an object returned by map_colors()", call. = FALSE)
  }
  n <- validate_nonneg_int(n, "n", 2L)
  n_ticks <- validate_nonneg_int(n_ticks, "n_ticks", 2L)
  cols <- grDevices::hcl.colors(n, palette = palette, alpha = alpha, rev = rev)
  ticks <- seq(limits[1L], limits[2L], length.out = n_ticks)
  tick_labels <- formatC(ticks, digits = digits, format = "fg")

  if (isTRUE(horizontal)) {
    op <- graphics::par(mar = c(2.6, 1.0, if (is.null(title)) 0.8 else 2.2, 1.0))
    on.exit(graphics::par(op), add = TRUE)
    graphics::plot.new()
    graphics::plot.window(xlim = limits, ylim = c(0, 1), xaxs = "i", yaxs = "i")
    xs <- seq(limits[1L], limits[2L], length.out = n + 1L)
    if (limits[1L] == limits[2L]) {
      graphics::rect(limits[1L] - 0.5, 0, limits[2L] + 0.5, 1, col = cols[ceiling(n / 2)], border = NA)
    } else {
      graphics::rect(xs[-length(xs)], 0, xs[-1L], 1, col = cols, border = NA)
    }
    graphics::box()
    graphics::axis(1L, at = ticks, labels = tick_labels, cex.axis = cex_axis, las = 1L, ...)
    if (!is.null(title)) graphics::title(main = title, cex.main = 0.95)
  } else {
    op <- graphics::par(mar = c(1.2, 0.8, if (is.null(title)) 1.0 else 2.2, 3.6))
    on.exit(graphics::par(op), add = TRUE)
    graphics::plot.new()
    graphics::plot.window(xlim = c(0, 1), ylim = limits, xaxs = "i", yaxs = "i")
    ys <- seq(limits[1L], limits[2L], length.out = n + 1L)
    if (limits[1L] == limits[2L]) {
      graphics::rect(0, limits[1L] - 0.5, 1, limits[2L] + 0.5, col = cols[ceiling(n / 2)], border = NA)
    } else {
      graphics::rect(0, ys[-length(ys)], 1, ys[-1L], col = cols, border = NA)
    }
    graphics::box()
    graphics::axis(4L, at = ticks, labels = tick_labels, cex.axis = cex_axis, las = 1L, ...)
    if (!is.null(title)) graphics::title(main = title, cex.main = 0.95)
  }
  invisible(limits)
}

resolve_plot_format <- function(output, fmt = "auto") {
  fmt <- tolower(as.character(fmt)[1L])
  if (fmt == "auto") {
    ext <- tolower(tools::file_ext(output))
    if (!nzchar(ext)) {
      stop("Cannot infer format from output filename without extension; specify fmt explicitly",
           call. = FALSE)
    }
    fmt <- switch(ext,
      pdf = "pdf",
      svg = "svg",
      png = "png",
      ps = "ps",
      eps = "eps",
      jpg = "jpeg",
      jpeg = "jpeg",
      tif = "tiff",
      tiff = "tiff",
      stop(sprintf("Unsupported output file extension '.%s'", ext), call. = FALSE)
    )
  }
  allowed <- c("pdf", "svg", "png", "ps", "eps", "jpeg", "tiff")
  if (!(fmt %in% allowed)) {
    stop(sprintf("fmt must be one of: %s", paste(allowed, collapse = ", ")), call. = FALSE)
  }
  fmt
}

#' @noRd
#' @export
save_graph_plot <- function(expr, output, width = 7, height = 5, dpi = 150,
                            fmt = c("auto", "pdf", "svg", "png", "ps", "eps", "jpeg", "tiff"),
                            bg = "white") {
  if (!is.character(output) || length(output) != 1L || !nzchar(output)) {
    stop("output must be a non-empty file path", call. = FALSE)
  }
  width <- validate_pos_scalar(width, "width")
  height <- validate_pos_scalar(height, "height")
  dpi <- validate_pos_scalar(dpi, "dpi")
  fmt <- resolve_plot_format(output, match.arg(fmt))
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)

  if (fmt == "pdf") {
    grDevices::pdf(file = output, width = width, height = height, bg = bg)
  } else if (fmt == "svg") {
    if (requireNamespace("svglite", quietly = TRUE)) {
      svglite::svglite(file = output, width = width, height = height, bg = bg)
    } else {
      ok_svg <- tryCatch({
        grDevices::svg(filename = output, width = width, height = height, bg = bg)
        TRUE
      }, error = function(e) FALSE, warning = function(w) {
        if (grepl("cairo", conditionMessage(w), ignore.case = TRUE)) FALSE else {
          warning(w)
          TRUE
        }
      })
      if (!isTRUE(ok_svg)) {
        stop("SVG device unavailable (Cairo/X11 library could not be loaded); install 'svglite' or export as 'pdf'/'eps'",
             call. = FALSE)
      }
    }
  } else if (fmt == "ps") {
    grDevices::postscript(file = output, width = width, height = height,
                          horizontal = FALSE, onefile = TRUE, paper = "special", bg = bg)
  } else if (fmt == "eps") {
    grDevices::postscript(file = output, width = width, height = height,
                          horizontal = FALSE, onefile = FALSE, paper = "special", bg = bg)
  } else if (fmt == "png") {
    grDevices::png(filename = output, width = width, height = height,
                   units = "in", res = dpi, bg = bg)
  } else if (fmt == "jpeg") {
    grDevices::jpeg(filename = output, width = width, height = height,
                    units = "in", res = dpi, bg = bg, quality = 92)
  } else if (fmt == "tiff") {
    grDevices::tiff(filename = output, width = width, height = height,
                    units = "in", res = dpi, bg = bg)
  }
  on.exit(grDevices::dev.off(), add = TRUE)

  val <- expr
  if (inherits(val, "ggplot")) {
    print(val)
  } else if (inherits(val, "grob")) {
    grid::grid.newpage()
    grid::grid.draw(val)
  }
  invisible(output)
}

quadratic_curve_points <- function(x0, y0, x1, y1, curvature = 0.15, n = 16L) {
  mx <- 0.5 * (x0 + x1)
  my <- 0.5 * (y0 + y1)
  dx <- x1 - x0
  dy <- y1 - y0
  cx <- mx - curvature * dy
  cy <- my + curvature * dx
  t <- seq(0, 1, length.out = n)
  mt <- 1 - t
  cbind(
    x = mt * mt * x0 + 2 * mt * t * cx + t * t * x1,
    y = mt * mt * y0 + 2 * mt * t * cy + t * t * y1
  )
}

#' @noRd
#' @export
graph_to_data_frames <- function(graph, pos = NULL, vertex_data = NULL, edge_data = NULL) {
  graph_pointer(graph)
  n <- graph_num_vertices(graph)
  if (is.null(pos)) {
    pos <- sfdp_layout(graph)
  } else {
    pos <- validate_pos_matrix(pos, n, "pos")
  }
  vdf <- data.frame(
    id = seq_len(n),
    x = if (n > 0L) pos[, 1L] else numeric(0),
    y = if (n > 0L) pos[, 2L] else numeric(0),
    degree = if (n > 0L) as.numeric(graph_degree(graph)) else numeric(0),
    stringsAsFactors = FALSE
  )
  if (!is.null(vertex_data)) {
    vd <- as.data.frame(vertex_data, stringsAsFactors = FALSE)
    if (nrow(vd) != n) {
      stop("vertex_data must have nrow equal to graph_num_vertices(graph)", call. = FALSE)
    }
    for (nm in names(vd)) vdf[[nm]] <- vd[[nm]]
  }

  edges <- graph_edges(graph)
  m <- nrow(edges)
  edf <- data.frame(
    from = if (m > 0L) edges[, 1L] else integer(0),
    to = if (m > 0L) edges[, 2L] else integer(0),
    x = if (m > 0L) pos[edges[, 1L], 1L] else numeric(0),
    y = if (m > 0L) pos[edges[, 1L], 2L] else numeric(0),
    xend = if (m > 0L) pos[edges[, 2L], 1L] else numeric(0),
    yend = if (m > 0L) pos[edges[, 2L], 2L] else numeric(0),
    stringsAsFactors = FALSE
  )
  eprops <- graph_edge_properties(graph)
  for (nm in eprops) {
    edf[[nm]] <- graph_edge_property(graph, nm)
  }
  if (!is.null(edge_data)) {
    ed <- as.data.frame(edge_data, stringsAsFactors = FALSE)
    if (nrow(ed) != m) {
      stop("edge_data must have nrow equal to graph_num_edges(graph)", call. = FALSE)
    }
    for (nm in names(ed)) edf[[nm]] <- ed[[nm]]
  }
  list(vertices = vdf, edges = edf)
}

#' @noRd
#' @export
graph_grob <- function(graph, pos = NULL,
                       vertex_color = "#4C72B0",
                       vertex_border = "#1F2937",
                       vertex_size = 1.0,
                       vertex_halo = FALSE,
                       vertex_halo_color = grDevices::adjustcolor("#E63946", 0.4),
                       vertex_halo_size = 1.8,
                       vertex_visible = TRUE,
                       vertex_label = NULL,
                       vertex_label_cex = 0.7,
                       edge_color = grDevices::adjustcolor("#4B5563", 0.45),
                       edge_lwd = 1.0,
                       edge_curve = 0,
                       edge_visible = TRUE,
                       vorder = NULL,
                       eorder = NULL,
                       arrows = NULL,
                       xlim = NULL,
                       ylim = NULL,
                       pad = 0.08,
                       bg = NULL,
                       axes = FALSE,
                       main = NULL) {
  graph_pointer(graph)
  n <- graph_num_vertices(graph)
  if (is.null(pos)) {
    pos <- sfdp_layout(graph)
  } else {
    pos <- validate_pos_matrix(pos, n, "pos")
  }
  if (is.null(arrows)) arrows <- isTRUE(graph$directed)
  lims <- compute_plot_limits(pos, pad = pad, xlim = xlim, ylim = ylim)

  vp <- grid::viewport(
    x = grid::unit(0.5, "npc"),
    y = grid::unit(0.5, "npc"),
    width = grid::unit(if (isTRUE(axes)) 0.84 else 0.94, "npc"),
    height = grid::unit(if (!is.null(main)) 0.84 else (if (isTRUE(axes)) 0.84 else 0.94), "npc"),
    xscale = lims$xlim,
    yscale = lims$ylim
  )

  children <- list()
  if (!is.null(bg)) {
    children[[length(children) + 1L]] <- grid::rectGrob(
      gp = grid::gpar(fill = bg, col = NA), vp = vp
    )
  }
  if (isTRUE(axes)) {
    children[[length(children) + 1L]] <- grid::xaxisGrob(vp = vp, gp = grid::gpar(cex = 0.75))
    children[[length(children) + 1L]] <- grid::yaxisGrob(vp = vp, gp = grid::gpar(cex = 0.75))
    children[[length(children) + 1L]] <- grid::rectGrob(
      gp = grid::gpar(fill = NA, col = "#64748B", lwd = 0.8), vp = vp
    )
  }

  v_vis <- rep_len(as.logical(vertex_visible), n)
  v_vis[is.na(v_vis)] <- FALSE
  edges <- graph_edges(graph)
  m <- nrow(edges)
  if (m > 0L) {
    e_vis <- rep_len(as.logical(edge_visible), m)
    e_vis[is.na(e_vis)] <- FALSE
    u <- edges[, 1L]
    v <- edges[, 2L]
    keep <- which(e_vis & v_vis[u] & v_vis[v] & (u != v))
    if (length(keep) > 0L) {
      if (!is.null(eorder)) {
        keep <- keep[order(eorder[keep])]
      }
      e_col <- rep_len(edge_color, m)[keep]
      e_lwd <- rep_len(edge_lwd, m)[keep]
      e_crv <- rep_len(as.numeric(edge_curve), m)[keep]
      arr_spec <- if (isTRUE(arrows)) {
        grid::arrow(angle = 20, length = grid::unit(0.08, "inches"), type = "closed")
      } else {
        NULL
      }
      if (all(e_crv == 0)) {
        children[[length(children) + 1L]] <- grid::segmentsGrob(
          x0 = grid::unit(pos[u[keep], 1L], "native"),
          y0 = grid::unit(pos[u[keep], 2L], "native"),
          x1 = grid::unit(pos[v[keep], 1L], "native"),
          y1 = grid::unit(pos[v[keep], 2L], "native"),
          arrow = arr_spec,
          gp = grid::gpar(col = e_col, fill = e_col, lwd = e_lwd),
          vp = vp
        )
      } else {
        for (idx in seq_along(keep)) {
          k <- keep[idx]
          pts <- quadratic_curve_points(pos[u[k], 1L], pos[u[k], 2L],
                                        pos[v[k], 1L], pos[v[k], 2L],
                                        curvature = e_crv[idx], n = 14L)
          children[[length(children) + 1L]] <- grid::linesGrob(
            x = grid::unit(pts[, 1L], "native"),
            y = grid::unit(pts[, 2L], "native"),
            arrow = arr_spec,
            gp = grid::gpar(col = e_col[idx], fill = e_col[idx], lwd = e_lwd[idx]),
            vp = vp
          )
        }
      }
    }
  }

  if (n > 0L) {
    v_col <- rep_len(vertex_color, n)
    v_brd <- rep_len(vertex_border, n)
    v_siz <- rep_len(vertex_size, n)
    v_hal <- rep_len(as.logical(vertex_halo), n)
    v_hal[is.na(v_hal)] <- FALSE
    v_hcol <- rep_len(vertex_halo_color, n)
    v_hsiz <- rep_len(vertex_halo_size, n)

    vis_idx <- which(v_vis)
    if (!is.null(vorder) && length(vis_idx) > 0L) {
      vis_idx <- vis_idx[order(vorder[vis_idx])]
    }
    hal_idx <- vis_idx[v_hal[vis_idx]]
    if (length(hal_idx) > 0L) {
      children[[length(children) + 1L]] <- grid::pointsGrob(
        x = grid::unit(pos[hal_idx, 1L], "native"),
        y = grid::unit(pos[hal_idx, 2L], "native"),
        pch = 19L,
        size = grid::unit(v_hsiz[hal_idx] * 0.45, "char"),
        gp = grid::gpar(col = v_hcol[hal_idx]),
        vp = vp
      )
    }
    if (length(vis_idx) > 0L) {
      children[[length(children) + 1L]] <- grid::pointsGrob(
        x = grid::unit(pos[vis_idx, 1L], "native"),
        y = grid::unit(pos[vis_idx, 2L], "native"),
        pch = 21L,
        size = grid::unit(v_siz[vis_idx] * 0.45, "char"),
        gp = grid::gpar(fill = v_col[vis_idx], col = v_brd[vis_idx], lwd = 0.8),
        vp = vp
      )
      if (!is.null(vertex_label)) {
        v_lab <- rep_len(as.character(vertex_label), n)
        children[[length(children) + 1L]] <- grid::textGrob(
          label = v_lab[vis_idx],
          x = grid::unit(pos[vis_idx, 1L], "native"),
          y = grid::unit(pos[vis_idx, 2L], "native"),
          gp = grid::gpar(cex = vertex_label_cex),
          vp = vp
        )
      }
    }
  }

  if (!is.null(main)) {
    children[[length(children) + 1L]] <- grid::textGrob(
      label = as.character(main)[1L],
      x = grid::unit(0.5, "npc"),
      y = grid::unit(0.96, "npc"),
      gp = grid::gpar(fontface = "bold", cex = 1.05)
    )
  }

  grid::gTree(children = do.call(grid::gList, children), cl = "rgraphtool_grob")
}

#' @noRd
#' @export
ggplot_graph <- function(graph, pos = NULL,
                         vertex_color = "#4C72B0",
                         vertex_size = 2.6,
                         vertex_label = NULL,
                         edge_color = "#64748B",
                         edge_alpha = 0.5,
                         edge_lwd = 0.5,
                         edge_curve = 0,
                         arrows = NULL,
                         palette = "Viridis",
                         color_label = NULL,
                         size_label = NULL,
                         axes = FALSE,
                         main = NULL,
                         subtitle = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for ggplot_graph(); install it or use plot_graph() / graph_grob()",
         call. = FALSE)
  }
  graph_pointer(graph)
  n <- graph_num_vertices(graph)
  dfs <- graph_to_data_frames(graph, pos = pos)
  vdf <- dfs$vertices
  edf <- dfs$edges[dfs$edges$from != dfs$edges$to, , drop = FALSE]
  if (is.null(arrows)) arrows <- isTRUE(graph$directed)

  p <- ggplot2::ggplot()
  x <- y <- xend <- yend <- .ecol <- .vfill <- .vsize <- .vlabel <- NULL
  arr_spec <- if (isTRUE(arrows)) {
    grid::arrow(angle = 20, length = grid::unit(0.08, "inches"), type = "closed")
  } else {
    NULL
  }

  if (nrow(edf) > 0L) {
    if (is.numeric(edge_color) && length(edge_color) == nrow(dfs$edges)) {
      edf$.ecol <- edge_color[dfs$edges$from != dfs$edges$to]
      if (edge_curve != 0) {
        p <- p + ggplot2::geom_curve(
          data = edf,
          mapping = ggplot2::aes(x = x, y = y, xend = xend, yend = yend, colour = .ecol),
          curvature = edge_curve, alpha = edge_alpha, linewidth = edge_lwd, arrow = arr_spec
        )
      } else {
        p <- p + ggplot2::geom_segment(
          data = edf,
          mapping = ggplot2::aes(x = x, y = y, xend = xend, yend = yend, colour = .ecol),
          alpha = edge_alpha, linewidth = edge_lwd, arrow = arr_spec
        )
      }
      p <- p + ggplot2::scale_colour_gradientn(
        colours = grDevices::hcl.colors(256L, palette = palette),
        name = if (is.null(color_label)) "Edge value" else color_label
      )
    } else {
      if (edge_curve != 0) {
        p <- p + ggplot2::geom_curve(
          data = edf,
          mapping = ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
          colour = edge_color[1L], alpha = edge_alpha, linewidth = edge_lwd,
          curvature = edge_curve, arrow = arr_spec
        )
      } else {
        p <- p + ggplot2::geom_segment(
          data = edf,
          mapping = ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
          colour = edge_color[1L], alpha = edge_alpha, linewidth = edge_lwd, arrow = arr_spec
        )
      }
    }
  }

  if (nrow(vdf) > 0L) {
    vcol_mapped <- (length(vertex_color) == n && (is.numeric(vertex_color) || is.factor(vertex_color)))
    vsiz_mapped <- (length(vertex_size) == n && is.numeric(vertex_size))

    if (vcol_mapped) vdf$.vfill <- vertex_color
    if (vsiz_mapped) vdf$.vsize <- vertex_size

    if (vcol_mapped && vsiz_mapped) {
      p <- p + ggplot2::geom_point(
        data = vdf,
        mapping = ggplot2::aes(x = x, y = y, fill = .vfill, size = .vsize),
        shape = 21L, colour = "#1F2937", stroke = 0.4
      )
    } else if (vcol_mapped) {
      p <- p + ggplot2::geom_point(
        data = vdf,
        mapping = ggplot2::aes(x = x, y = y, fill = .vfill),
        size = vertex_size[1L], shape = 21L, colour = "#1F2937", stroke = 0.4
      )
    } else if (vsiz_mapped) {
      p <- p + ggplot2::geom_point(
        data = vdf,
        mapping = ggplot2::aes(x = x, y = y, size = .vsize),
        fill = vertex_color[1L], shape = 21L, colour = "#1F2937", stroke = 0.4
      )
    } else {
      vcol_vec <- rep_len(as.character(vertex_color), n)
      vsiz_vec <- rep_len(as.numeric(vertex_size), n)
      p <- p + ggplot2::geom_point(
        data = vdf,
        mapping = ggplot2::aes(x = x, y = y),
        fill = vcol_vec, size = vsiz_vec, shape = 21L, colour = "#1F2937", stroke = 0.4
      )
    }

    if (vcol_mapped) {
      if (is.numeric(vertex_color)) {
        p <- p + ggplot2::scale_fill_gradientn(
          colours = grDevices::hcl.colors(256L, palette = palette),
          name = if (is.null(color_label)) "Value" else color_label
        )
      } else {
        p <- p + ggplot2::scale_fill_discrete(
          name = if (is.null(color_label)) "Group" else color_label
        )
      }
    }
    if (vsiz_mapped) {
      p <- p + ggplot2::scale_size_continuous(
        range = c(1.8, 6.0),
        name = if (is.null(size_label)) "Size" else size_label
      )
    }
    if (!is.null(vertex_label)) {
      vdf$.vlabel <- rep_len(as.character(vertex_label), n)
      p <- p + ggplot2::geom_text(
        data = vdf,
        mapping = ggplot2::aes(x = x, y = y, label = .vlabel),
        size = 2.6, vjust = -0.9, check_overlap = TRUE
      )
    }
  }

  p <- p + ggplot2::coord_equal()
  if (isTRUE(axes)) {
    p <- p + ggplot2::theme_minimal(base_size = 11)
  } else {
    p <- p + ggplot2::theme_void(base_size = 11) +
      ggplot2::theme(plot.margin = ggplot2::margin(6, 6, 6, 6))
  }
  if (!is.null(main) || !is.null(subtitle)) {
    p <- p + ggplot2::labs(title = main, subtitle = subtitle)
  }
  p
}
