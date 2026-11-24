# Geographic map projections, great-circle geodesic edges, offline world
# coastlines, and geographic network visualization (Basemap/Cartopy equivalent).

#' @noRd
#' @export
project_coords <- function(lon, lat,
                           projection = c("robinson", "mercator", "equirectangular",
                                          "mollweide", "orthographic"),
                           lon0 = 0, lat0 = 0) {
  projection <- match.arg(projection)
  if (!is.numeric(lon) || !is.numeric(lat) || length(lon) != length(lat)) {
    stop("lon and lat must be numeric vectors of equal length", call. = FALSE)
  }
  if (length(lon) == 0L) {
    return(data.frame(x = numeric(0), y = numeric(0), visible = logical(0)))
  }
  # Wrap relative longitude to [-180, 180]
  dlon <- ((lon - lon0 + 180) %% 360) - 180
  rad <- pi / 180
  lam <- dlon * rad
  phi <- pmax(-90, pmin(90, lat)) * rad
  visible <- rep(TRUE, length(lon))

  if (projection == "equirectangular") {
    x <- lam * cos(lat0 * rad)
    y <- phi
  } else if (projection == "mercator") {
    phi_c <- pmax(-85, pmin(85, lat)) * rad
    x <- lam
    y <- log(tan(pi / 4 + phi_c / 2))
  } else if (projection == "robinson") {
    # Standard Robinson (1974) table at 5-degree latitude intervals (0..90 deg)
    lat_tbl <- seq(0, 90, by = 5)
    plen <- c(1.0000, 0.9986, 0.9954, 0.9900, 0.9822, 0.9730, 0.9600,
              0.9427, 0.9216, 0.8962, 0.8679, 0.8350, 0.7986, 0.7597,
              0.7186, 0.6732, 0.6213, 0.5722, 0.5322)
    pdfe <- c(0.0000, 0.0620, 0.1240, 0.1860, 0.2480, 0.3100, 0.3720,
              0.4340, 0.4958, 0.5571, 0.6176, 0.6769, 0.7346, 0.7903,
              0.8435, 0.8936, 0.9394, 0.9761, 1.0000)
    alat <- abs(pmax(-90, pmin(90, lat)))
    px <- stats::approx(lat_tbl, plen, xout = alat, rule = 2)$y
    py <- stats::approx(lat_tbl, pdfe, xout = alat, rule = 2)$y
    x <- 0.8487 * px * lam
    y <- 1.3523 * sign(lat) * py
  } else if (projection == "mollweide") {
    theta <- phi
    target <- pi * sin(phi)
    for (iter in seq_len(8L)) {
      f <- 2 * theta + sin(2 * theta) - target
      df <- 2 + 2 * cos(2 * theta)
      step <- ifelse(abs(df) > 1e-12, f / df, 0)
      theta <- theta - step
    }
    x <- (2 * sqrt(2) / pi) * lam * cos(theta)
    y <- sqrt(2) * sin(theta)
  } else if (projection == "orthographic") {
    phi0 <- lat0 * rad
    cos_c <- sin(phi0) * sin(phi) + cos(phi0) * cos(phi) * cos(lam)
    x <- cos(phi) * sin(lam)
    y <- cos(phi0) * sin(phi) - sin(phi0) * cos(phi) * cos(lam)
    visible <- (cos_c >= -1e-9)
  }

  visible[is.na(lon) | is.na(lat)] <- FALSE
  data.frame(x = x, y = y, visible = visible)
}

#' @noRd
#' @export
great_circle_arc <- function(lon1, lat1, lon2, lat2, n = 24L) {
  n <- validate_nonneg_int(n, "n", 2L)
  rad <- pi / 180
  p1 <- c(cos(lat1 * rad) * cos(lon1 * rad),
          cos(lat1 * rad) * sin(lon1 * rad),
          sin(lat1 * rad))
  p2 <- c(cos(lat2 * rad) * cos(lon2 * rad),
          cos(lat2 * rad) * sin(lon2 * rad),
          sin(lat2 * rad))
  dot <- max(-1, min(1, sum(p1 * p2)))
  omega <- acos(dot)
  t <- seq(0, 1, length.out = n)
  if (omega < 1e-6) {
    lons <- lon1 + t * (lon2 - lon1)
    lats <- lat1 + t * (lat2 - lat1)
  } else if (abs(omega - pi) < 1e-6) {
    # Antipodal fallback via small perturbation
    lons <- lon1 + t * (lon2 - lon1)
    lats <- lat1 + t * (lat2 - lat1)
  } else {
    s <- sin(omega)
    w1 <- sin((1 - t) * omega) / s
    w2 <- sin(t * omega) / s
    px <- w1 * p1[1L] + w2 * p2[1L]
    py <- w1 * p1[2L] + w2 * p2[2L]
    pz <- w1 * p1[3L] + w2 * p2[3L]
    lons <- atan2(py, px) / rad
    lats <- asin(pmax(-1, pmin(1, pz))) / rad
  }
  jumps <- c(FALSE, abs(diff(lons)) > 180)
  seg <- cumsum(jumps) + 1L
  data.frame(lon = lons, lat = lats, segment = as.integer(seg))
}

#' @noRd
#' @export
world_coastlines <- function() {
  polys <- list(
    north_america = rbind(
      c(-166, 60), c(-162, 66), c(-156, 71), c(-138, 69), c(-124, 70),
      c(-104, 68), c(-88, 69), c(-78, 63), c(-64, 60), c(-55, 52),
      c(-60, 46), c(-67, 45), c(-70, 42), c(-75, 39), c(-76, 35),
      c(-81, 31), c(-80, 25), c(-83, 29), c(-89, 30), c(-97, 26),
      c(-96, 19), c(-88, 21), c(-83, 15), c(-83, 9), c(-78, 8),
      c(-84, 10), c(-92, 15), c(-105, 20), c(-110, 24), c(-114, 31),
      c(-117, 32), c(-121, 35), c(-124, 41), c(-124, 48), c(-131, 54),
      c(-140, 60), c(-152, 60), c(-166, 60)
    ),
    greenland = rbind(
      c(-44, 60), c(-52, 65), c(-54, 72), c(-66, 77), c(-48, 82),
      c(-22, 82), c(-18, 75), c(-22, 70), c(-36, 65), c(-44, 60)
    ),
    south_america = rbind(
      c(-78, 9), c(-72, 11), c(-62, 10), c(-50, 1), c(-35, -6),
      c(-38, -14), c(-42, -23), c(-49, -28), c(-54, -35), c(-62, -40),
      c(-66, -47), c(-69, -54), c(-74, -52), c(-73, -42), c(-71, -30),
      c(-70, -18), c(-77, -12), c(-81, -5), c(-78, 2), c(-78, 9)
    ),
    europe_asia = rbind(
      c(-9, 37), c(-9, 43), c(-2, 43), c(-4, 48), c(2, 51),
      c(8, 54), c(10, 58), c(5, 61), c(12, 66), c(20, 70),
      c(31, 70), c(44, 67), c(60, 69), c(72, 73), c(104, 77),
      c(135, 72), c(170, 69), c(179, 66), c(162, 60), c(156, 51),
      c(141, 53), c(135, 44), c(129, 35), c(126, 38), c(122, 37),
      c(121, 31), c(117, 24), c(109, 20), c(108, 12), c(104, 2),
      c(100, 13), c(97, 16), c(89, 22), c(80, 13), c(77, 8),
      c(72, 20), c(67, 25), c(56, 26), c(59, 22), c(50, 15),
      c(43, 13), c(35, 29), c(36, 36), c(27, 37), c(24, 38),
      c(20, 40), c(14, 45), c(16, 40), c(12, 42), c(8, 44),
      c(3, 42), c(-5, 36), c(-9, 37)
    ),
    uk = rbind(
      c(-6, 50), c(-5, 55), c(-3, 58), c(1, 53), c(1, 51), c(-6, 50)
    ),
    africa = rbind(
      c(-6, 35), c(10, 37), c(15, 32), c(25, 32), c(32, 31),
      c(35, 27), c(43, 12), c(51, 11), c(45, 2), c(40, -5),
      c(40, -15), c(35, -24), c(28, -33), c(18, -34), c(15, -26),
      c(12, -16), c(12, -6), c(9, 4), c(2, 6), c(-8, 5),
      c(-15, 11), c(-17, 15), c(-16, 21), c(-10, 30), c(-6, 35)
    ),
    madagascar = rbind(
      c(49, -12), c(50, -17), c(47, -25), c(44, -24), c(44, -17), c(49, -12)
    ),
    japan = rbind(
      c(130, 32), c(135, 34), c(140, 36), c(142, 41), c(145, 44),
      c(140, 42), c(136, 37), c(131, 34), c(130, 32)
    ),
    australia = rbind(
      c(114, -22), c(122, -17), c(130, -12), c(136, -12), c(142, -11),
      c(146, -19), c(153, -27), c(150, -37), c(141, -38), c(135, -34),
      c(126, -33), c(115, -34), c(113, -26), c(114, -22)
    ),
    new_zealand = rbind(
      c(173, -35), c(178, -38), c(175, -41), c(171, -44), c(167, -46),
      c(170, -42), c(174, -39), c(173, -35)
    ),
    antarctica = rbind(
      c(-180, -72), c(-120, -74), c(-65, -65), c(-55, -72), c(-20, -73),
      c(20, -70), c(70, -68), c(120, -67), c(165, -72), c(180, -72),
      c(180, -85), c(-180, -85), c(-180, -72)
    )
  )
  out <- vector("list", length(polys))
  for (i in seq_along(polys)) {
    mat <- polys[[i]]
    out[[i]] <- data.frame(
      lon = mat[, 1L],
      lat = mat[, 2L],
      group = rep(as.integer(i), nrow(mat)),
      region = rep(names(polys)[i], nrow(mat)),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

#' @noRd
#' @export
plot_graph_map <- function(graph, lon, lat,
                           projection = c("robinson", "mercator", "equirectangular",
                                          "mollweide", "orthographic"),
                           lon0 = 0, lat0 = 0,
                           xlim = NULL, ylim = NULL,
                           geodesic = TRUE, n_arc = 24L,
                           coastlines = TRUE, polygons = NULL,
                           graticule = TRUE,
                           land_col = "#ECEAE4",
                           border_col = "#94A3B8",
                           ocean_col = "#E0F2FE",
                           graticule_col = grDevices::adjustcolor("#64748B", 0.30),
                           vertex_color = "#D62828",
                           vertex_border = "#1F2937",
                           vertex_size = 1.2,
                           vertex_label = NULL,
                           vertex_label_cex = 0.7,
                           vertex_label_pos = 3L,
                           edge_color = grDevices::adjustcolor("#1D4ED8", 0.55),
                           edge_lwd = 1.1,
                           eorder = NULL,
                           vorder = NULL,
                           arrows = NULL,
                           arrow_length = 0.06,
                           main = NULL,
                           output = NULL,
                           width = 9,
                           height = 5,
                           dpi = 150,
                           ...) {
  graph_pointer(graph)
  projection <- match.arg(projection)
  n <- graph_num_vertices(graph)
  if (!is.numeric(lon) || !is.numeric(lat) || length(lon) != n || length(lat) != n ||
      anyNA(lon) || anyNA(lat) || any(!is.finite(lon)) || any(!is.finite(lat))) {
    stop("lon and lat must be finite numeric vectors of length graph_num_vertices(graph)",
         call. = FALSE)
  }
  if (is.null(arrows)) arrows <- isTRUE(graph$directed)

  proj_v <- project_coords(lon, lat, projection = projection, lon0 = lon0, lat0 = lat0)

  render_map <- function() {
    # Determine projected plot limits
    if (!is.null(xlim) && !is.null(ylim)) {
      box_lons <- c(seq(xlim[1L], xlim[2L], length.out = 25L),
                    rep(xlim[2L], 25L),
                    seq(xlim[2L], xlim[1L], length.out = 25L),
                    rep(xlim[1L], 25L))
      box_lats <- c(rep(ylim[1L], 25L),
                    seq(ylim[1L], ylim[2L], length.out = 25L),
                    rep(ylim[2L], 25L),
                    seq(ylim[2L], ylim[1L], length.out = 25L))
      bp <- project_coords(box_lons, box_lats, projection = projection, lon0 = lon0, lat0 = lat0)
      pxlim <- range(bp$x[bp$visible], finite = TRUE)
      pylim <- range(bp$y[bp$visible], finite = TRUE)
    } else if (projection == "orthographic") {
      pxlim <- c(-1.05, 1.05)
      pylim <- c(-1.05, 1.05)
    } else {
      lat_bnd <- if (projection == "mercator") 80 else 89.9
      top_b <- project_coords(seq(-180, 180, length.out = 60L), rep(lat_bnd, 60L),
                              projection = projection, lon0 = lon0, lat0 = lat0)
      bot_b <- project_coords(seq(-180, 180, length.out = 60L), rep(-lat_bnd, 60L),
                              projection = projection, lon0 = lon0, lat0 = lat0)
      l_b <- project_coords(rep(lon0 - 179.9, 60L), seq(-lat_bnd, lat_bnd, length.out = 60L),
                            projection = projection, lon0 = lon0, lat0 = lat0)
      r_b <- project_coords(rep(lon0 + 179.9, 60L), seq(-lat_bnd, lat_bnd, length.out = 60L),
                            projection = projection, lon0 = lon0, lat0 = lat0)
      pxlim <- range(c(top_b$x, bot_b$x, l_b$x, r_b$x), finite = TRUE)
      pylim <- range(c(top_b$y, bot_b$y, l_b$y, r_b$y), finite = TRUE)
    }

    op <- graphics::par(bg = "white", mar = c(1.0, 1.0, if (is.null(main)) 1.0 else 2.4, 1.0))
    on.exit(graphics::par(op), add = TRUE)
    graphics::plot.new()
    graphics::plot.window(xlim = pxlim, ylim = pylim, asp = 1)

    # Ocean backdrop
    if (!is.null(ocean_col)) {
      if (projection == "orthographic") {
        ang <- seq(0, 2 * pi, length.out = 180L)
        graphics::polygon(cos(ang), sin(ang), col = ocean_col, border = border_col, lwd = 1.0)
      } else if (is.null(xlim) && is.null(ylim) && projection %in% c("robinson", "mollweide")) {
        lat_seq <- seq(-89.9, 89.9, length.out = 90L)
        left_p <- project_coords(rep(lon0 - 179.9, 90L), lat_seq, projection, lon0, lat0)
        right_p <- project_coords(rep(lon0 + 179.9, 90L), rev(lat_seq), projection, lon0, lat0)
        graphics::polygon(c(left_p$x, right_p$x), c(left_p$y, right_p$y),
                          col = ocean_col, border = border_col, lwd = 0.9)
      } else {
        graphics::rect(pxlim[1L], pylim[1L], pxlim[2L], pylim[2L],
                       col = ocean_col, border = border_col, lwd = 0.9)
      }
    }

    # Graticules (lat/lon grid)
    if (isTRUE(graticule)) {
      lat_max <- if (projection == "mercator") 80 else 85
      for (la in seq(-60, 60, by = 30)) {
        ln_seq <- seq(lon0 - 179.5, lon0 + 179.5, length.out = 120L)
        gp <- project_coords(ln_seq, rep(la, length(ln_seq)), projection, lon0, lat0)
        gx <- ifelse(gp$visible, gp$x, NA_real_)
        gy <- ifelse(gp$visible, gp$y, NA_real_)
        graphics::lines(gx, gy, col = graticule_col, lty = 3L, lwd = 0.7)
      }
      for (ln in seq(-150, 180, by = 30)) {
        la_seq <- seq(-lat_max, lat_max, length.out = 90L)
        gp <- project_coords(rep(ln, length(la_seq)), la_seq, projection, lon0, lat0)
        gx <- ifelse(gp$visible, gp$x, NA_real_)
        gy <- ifelse(gp$visible, gp$y, NA_real_)
        graphics::lines(gx, gy, col = graticule_col, lty = 3L, lwd = 0.7)
      }
    }

    # Land polygons / coastlines
    if (isTRUE(coastlines)) {
      poly_df <- if (is.null(polygons)) world_coastlines() else polygons
      if (!is.data.frame(poly_df) || !all(c("lon", "lat") %in% names(poly_df))) {
        stop("polygons must be a data.frame with columns 'lon' and 'lat' (and optional 'group')",
             call. = FALSE)
      }
      grp_col <- if ("group" %in% names(poly_df)) poly_df$group else rep(1L, nrow(poly_df))
      for (g_id in unique(grp_col)) {
        sub <- poly_df[grp_col == g_id, , drop = FALSE]
        pp <- project_coords(sub$lon, sub$lat, projection, lon0, lat0)
        # Skip polygons that wrap around the dateline in projected x or are completely hidden
        vis <- pp$visible
        if (sum(vis) < 3L) next
        if (projection == "orthographic") {
          px <- ifelse(vis, pp$x, NA_real_)
          py <- ifelse(vis, pp$y, NA_real_)
          if (all(vis)) {
            graphics::polygon(px, py, col = land_col, border = border_col, lwd = 0.7)
          } else {
            graphics::lines(px, py, col = border_col, lwd = 0.8)
          }
        } else {
          dl <- abs(diff(((sub$lon - lon0 + 180) %% 360) - 180))
          if (any(dl > 180, na.rm = TRUE)) {
            brk <- c(FALSE, dl > 180)
            seg_id <- cumsum(brk)
            for (sid in unique(seg_id)) {
              idx <- which(seg_id == sid)
              if (length(idx) >= 2L) {
                graphics::lines(pp$x[idx], pp$y[idx], col = border_col, lwd = 0.8)
              }
            }
          } else {
            graphics::polygon(pp$x, pp$y, col = land_col, border = border_col, lwd = 0.7)
          }
        }
      }
    }

    # Edges (geodesic great-circle arcs or projected straight lines)
    edges <- graph_edges(graph)
    m <- nrow(edges)
    if (m > 0L) {
      keep <- which(edges[, 1L] != edges[, 2L])
      if (length(keep) > 0L) {
        if (!is.null(eorder)) keep <- keep[order(eorder[keep])]
        e_col <- rep_len(edge_color, m)
        e_lwd <- rep_len(edge_lwd, m)
        for (k in keep) {
          u <- edges[k, 1L]
          v <- edges[k, 2L]
          if (isTRUE(geodesic)) {
            arc <- great_circle_arc(lon[u], lat[u], lon[v], lat[v], n = n_arc)
            for (sid in unique(arc$segment)) {
              sub <- arc[arc$segment == sid, , drop = FALSE]
              pa <- project_coords(sub$lon, sub$lat, projection, lon0, lat0)
              # Split if wrapping across lon0 +/- 180
              dl <- abs(diff(((sub$lon - lon0 + 180) %% 360) - 180))
              pa_x <- ifelse(pa$visible, pa$x, NA_real_)
              pa_y <- ifelse(pa$visible, pa$y, NA_real_)
              if (any(dl > 180)) {
                sub_seg <- cumsum(c(FALSE, dl > 180))
                for (ss in unique(sub_seg)) {
                  idx <- which(sub_seg == ss)
                  if (length(idx) >= 2L) {
                    graphics::lines(pa_x[idx], pa_y[idx], col = e_col[k], lwd = e_lwd[k])
                  }
                }
              } else {
                graphics::lines(pa_x, pa_y, col = e_col[k], lwd = e_lwd[k])
              }
            }
            if (isTRUE(arrows) && nrow(arc) >= 2L) {
              tail_sub <- arc[(nrow(arc) - 1L):nrow(arc), , drop = FALSE]
              pt <- project_coords(tail_sub$lon, tail_sub$lat, projection, lon0, lat0)
              if (all(pt$visible) && abs(pt$x[2L] - pt$x[1L]) + abs(pt$y[2L] - pt$y[1L]) > 1e-6) {
                graphics::arrows(pt$x[1L], pt$y[1L], pt$x[2L], pt$y[2L],
                                 col = e_col[k], lwd = e_lwd[k],
                                 length = arrow_length, angle = 20)
              }
            }
          } else {
            if (proj_v$visible[u] && proj_v$visible[v]) {
              if (isTRUE(arrows)) {
                graphics::arrows(proj_v$x[u], proj_v$y[u], proj_v$x[v], proj_v$y[v],
                                 col = e_col[k], lwd = e_lwd[k],
                                 length = arrow_length, angle = 20)
              } else {
                graphics::segments(proj_v$x[u], proj_v$y[u], proj_v$x[v], proj_v$y[v],
                                   col = e_col[k], lwd = e_lwd[k])
              }
            }
          }
        }
      }
    }

    # Vertices
    if (n > 0L) {
      vis_idx <- which(proj_v$visible)
      if (!is.null(vorder) && length(vis_idx) > 0L) {
        vis_idx <- vis_idx[order(vorder[vis_idx])]
      }
      if (length(vis_idx) > 0L) {
        v_col <- rep_len(vertex_color, n)[vis_idx]
        v_brd <- rep_len(vertex_border, n)[vis_idx]
        v_siz <- rep_len(vertex_size, n)[vis_idx]
        graphics::points(proj_v$x[vis_idx], proj_v$y[vis_idx],
                         pch = 21L, bg = v_col, col = v_brd, cex = v_siz, ...)
        if (!is.null(vertex_label)) {
          v_lab <- rep_len(as.character(vertex_label), n)[vis_idx]
          graphics::text(proj_v$x[vis_idx], proj_v$y[vis_idx],
                         labels = v_lab, cex = vertex_label_cex, pos = vertex_label_pos)
        }
      }
    }

    if (!is.null(main)) graphics::title(main = main)
  }

  if (!is.null(output)) {
    save_graph_plot(render_map(), output = output, width = width, height = height, dpi = dpi)
  } else {
    render_map()
  }

  invisible(list(projected = proj_v, projection = projection, lon0 = lon0, lat0 = lat0))
}
