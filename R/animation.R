# Layout, static plotting, animations, and interactive visualization
# R-native equivalents of graph-tool's animation & interactive cookbook demos.

validate_pos_matrix <- function(pos, n, arg = "pos") {
  if (!is.matrix(pos) || !is.numeric(pos) || ncol(pos) != 2L || nrow(pos) != n) {
    stop(sprintf("%s must be an %d x 2 numeric matrix", arg, n), call. = FALSE)
  }
  if (any(!is.finite(pos))) {
    stop(sprintf("%s must contain only finite coordinates", arg), call. = FALSE)
  }
  storage.mode(pos) <- "double"
  colnames(pos) <- c("x", "y")
  pos
}

validate_prob_scalar <- function(val, name) {
  if (!is.numeric(val) || length(val) != 1L || is.na(val) || !is.finite(val) ||
      val < 0 || val > 1) {
    stop(sprintf("%s must be a finite probability in [0, 1]", name), call. = FALSE)
  }
  as.double(val)
}

validate_pos_scalar <- function(val, name, allow_zero = FALSE) {
  if (!is.numeric(val) || length(val) != 1L || is.na(val) || !is.finite(val) ||
      (if (allow_zero) val < 0 else val <= 0)) {
    stop(sprintf("%s must be a %s finite number", name,
                 if (allow_zero) "non-negative" else "positive"), call. = FALSE)
  }
  as.double(val)
}

validate_nonneg_int <- function(val, name, min_val = 0L) {
  if (!is.numeric(val) || length(val) != 1L || is.na(val) || !is.finite(val) ||
      val < min_val || val != floor(val)) {
    stop(sprintf("%s must be an integer >= %d", name, min_val), call. = FALSE)
  }
  as.integer(val)
}

html_escape <- function(x) {
  x <- as.character(x)[1L]
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

json_num_vec <- function(x, digits = 5L) {
  if (length(x) == 0L) return("[]")
  fmt <- paste0("%.", digits, "g")
  paste0("[", paste(sprintf(fmt, as.numeric(x)), collapse = ","), "]")
}

json_int_vec <- function(x) {
  if (length(x) == 0L) return("[]")
  paste0("[", paste(as.integer(x), collapse = ","), "]")
}

json_pos_matrix <- function(pos) {
  n <- nrow(pos)
  if (n == 0L) return("[]")
  pairs <- sprintf("[%.5g,%.5g]", pos[, 1L], pos[, 2L])
  paste0("[", paste(pairs, collapse = ","), "]")
}

json_edge_matrix <- function(edges) {
  m <- nrow(edges)
  if (m == 0L) return("[]")
  pairs <- sprintf("[%d,%d]", as.integer(edges[, 1L]) - 1L, as.integer(edges[, 2L]) - 1L)
  paste0("[", paste(pairs, collapse = ","), "]")
}

sfdp_layout <- function(graph, pos = NULL, C = 0.2, K = 1.0, p = 2.0,
                        theta = 0.6, gamma = 1.0, mu = 0.0, mu_p = 1.0,
                        init_step = NULL, cooling_step = 0.95,
                        max_level = 10L, tol = 1e-4, max_iter = 200L,
                        simple = FALSE, groups = NULL, pin = NULL) {
  graph_pointer(graph)
  n <- graph_num_vertices(graph)
  if (!is.null(pos)) {
    pos <- validate_pos_matrix(pos, n, "pos")
  }
  C <- validate_pos_scalar(C, "C")
  K <- validate_pos_scalar(K, "K")
  p <- validate_pos_scalar(p, "p")
  theta <- validate_pos_scalar(theta, "theta", allow_zero = TRUE)
  if (!is.numeric(gamma) || length(gamma) != 1L || is.na(gamma) || !is.finite(gamma)) {
    stop("gamma must be a finite numeric scalar", call. = FALSE)
  }
  mu <- validate_pos_scalar(mu, "mu", allow_zero = TRUE)
  if (!is.numeric(mu_p) || length(mu_p) != 1L || is.na(mu_p) || !is.finite(mu_p)) {
    stop("mu_p must be a finite numeric scalar", call. = FALSE)
  }
  step_val <- if (is.null(init_step)) -1.0 else validate_pos_scalar(init_step, "init_step")
  cooling_step <- validate_pos_scalar(cooling_step, "cooling_step")
  if (cooling_step > 1.0) {
    stop("cooling_step must be in (0, 1]", call. = FALSE)
  }
  max_level <- validate_nonneg_int(max_level, "max_level", 0L)
  tol <- validate_pos_scalar(tol, "tol", allow_zero = TRUE)
  max_iter <- validate_nonneg_int(max_iter, "max_iter", 0L)
  if (!is.logical(simple) || length(simple) != 1L || is.na(simple)) {
    stop("simple must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.null(groups)) {
    if (length(groups) != n || anyNA(groups)) {
      stop("groups must have length equal to graph_num_vertices(graph) with no NA", call. = FALSE)
    }
    groups <- as.integer(factor(groups))
  }
  if (!is.null(pin)) {
    if (!is.logical(pin) || length(pin) != n || anyNA(pin)) {
      stop("pin must be a logical vector of length graph_num_vertices(graph) with no NA", call. = FALSE)
    }
  }
  cpp_sfdp_layout(
    pointer = graph$pointer,
    init_pos = pos,
    C = C,
    K = K,
    p = p,
    theta = theta,
    gamma = as.double(gamma),
    mu = mu,
    mu_p = as.double(mu_p),
    init_step = step_val,
    cooling_step = cooling_step,
    max_level = max_level,
    tol = tol,
    max_iter = max_iter,
    simple = simple,
    groups = groups,
    pin = pin
  )
}

graph_bfs <- function(graph, source, max_dist = Inf, directed = NULL) {
  graph_pointer(graph)
  n <- graph_num_vertices(graph)
  if (n == 0L) stop("graph has no vertices", call. = FALSE)
  source <- validate_nonneg_int(source, "source", 1L)
  if (source > n) {
    stop(sprintf("source must be between 1 and %d", n), call. = FALSE)
  }
  if (is.infinite(max_dist) && max_dist > 0) {
    md_int <- -1L
  } else {
    md_int <- validate_nonneg_int(max_dist, "max_dist", 0L)
  }
  if (is.null(directed)) {
    directed <- isTRUE(graph$directed)
  } else if (!is.logical(directed) || length(directed) != 1L || is.na(directed)) {
    stop("directed must be TRUE, FALSE, or NULL", call. = FALSE)
  }
  res <- cpp_graph_bfs(graph$pointer, source = source, max_dist = md_int, directed = directed)
  structure(res, class = "rgraphtool_bfs")
}

print.rgraphtool_bfs <- function(x, ...) {
  n <- length(x$distance)
  reached <- sum(is.finite(x$distance))
  max_d <- if (reached > 0L) max(x$distance[is.finite(x$distance)]) else NA_real_
  cat("rgraphtool BFS traversal\n")
  cat(sprintf("  source: %d (%s, max_dist = %s)\n",
              x$source,
              if (isTRUE(x$directed)) "directed" else "undirected",
              if (is.infinite(x$max_dist)) "Inf" else as.character(x$max_dist)))
  cat(sprintf("  reached: %d / %d vertices (max distance = %s, tree edges = %d)\n",
              reached, n, if (is.na(max_d)) "NA" else as.character(max_d),
              sum(x$tree_edges)))
  invisible(x)
}

compute_plot_limits <- function(pos, pad = 0.08, xlim = NULL, ylim = NULL) {
  if (nrow(pos) == 0L) {
    return(list(xlim = if (is.null(xlim)) c(-1, 1) else xlim,
                ylim = if (is.null(ylim)) c(-1, 1) else ylim))
  }
  xr <- range(pos[, 1L], finite = TRUE)
  yr <- range(pos[, 2L], finite = TRUE)
  dx <- diff(xr)
  dy <- diff(yr)
  span <- max(dx, dy)
  if (!is.finite(span) || span <= 0) span <- 1.0
  margin <- pad * span
  if (dx <= 0) xr <- xr + c(-0.5, 0.5)
  if (dy <= 0) yr <- yr + c(-0.5, 0.5)
  list(
    xlim = if (is.null(xlim)) c(xr[1L] - margin, xr[2L] + margin) else xlim,
    ylim = if (is.null(ylim)) c(yr[1L] - margin, yr[2L] + margin) else ylim
  )
}

plot_graph <- function(graph, pos = NULL,
                       vertex_color = "#4C72B0",
                       vertex_border = "#1F2937",
                       vertex_size = 1.2,
                       vertex_shape = 21L,
                       vertex_halo = FALSE,
                       vertex_halo_color = grDevices::adjustcolor("#E63946", 0.45),
                       vertex_halo_size = 2.2,
                       vertex_visible = TRUE,
                       vertex_label = NULL,
                       vertex_label_cex = 0.7,
                       edge_color = grDevices::adjustcolor("#4B5563", 0.45),
                       edge_lwd = 1.0,
                       edge_visible = TRUE,
                       arrows = NULL,
                       arrow_length = 0.08,
                       xlim = NULL,
                       ylim = NULL,
                       pad = 0.08,
                       bg = "white",
                       main = NULL,
                       add = FALSE,
                       ...) {
  graph_pointer(graph)
  n <- graph_num_vertices(graph)
  if (is.null(pos)) {
    pos <- sfdp_layout(graph)
  } else {
    pos <- validate_pos_matrix(pos, n, "pos")
  }
  if (is.null(arrows)) {
    arrows <- isTRUE(graph$directed)
  }
  edges <- graph_edges(graph)
  m <- nrow(edges)

  v_vis <- rep_len(as.logical(vertex_visible), n)
  v_vis[is.na(v_vis)] <- FALSE
  v_col <- rep_len(vertex_color, n)
  v_brd <- rep_len(vertex_border, n)
  v_siz <- rep_len(vertex_size, n)
  v_shp <- rep_len(vertex_shape, n)
  v_hal <- rep_len(as.logical(vertex_halo), n)
  v_hal[is.na(v_hal)] <- FALSE
  v_hcol <- rep_len(vertex_halo_color, n)
  v_hsiz <- rep_len(vertex_halo_size, n)

  if (!add) {
    lims <- compute_plot_limits(pos, pad = pad, xlim = xlim, ylim = ylim)
    op <- graphics::par(bg = bg, mar = c(1, 1, if (is.null(main)) 1 else 2.5, 1))
    on.exit(graphics::par(op), add = TRUE)
    graphics::plot.new()
    graphics::plot.window(xlim = lims$xlim, ylim = lims$ylim, asp = 1)
  }

  if (m > 0L) {
    e_vis <- rep_len(as.logical(edge_visible), m)
    e_vis[is.na(e_vis)] <- FALSE
    u <- edges[, 1L]
    v <- edges[, 2L]
    keep <- e_vis & v_vis[u] & v_vis[v] & (u != v)
    if (any(keep)) {
      e_col <- rep_len(edge_color, m)[keep]
      e_lwd <- rep_len(edge_lwd, m)[keep]
      x0 <- pos[u[keep], 1L]
      y0 <- pos[u[keep], 2L]
      x1 <- pos[v[keep], 1L]
      y1 <- pos[v[keep], 2L]
      if (isTRUE(arrows)) {
        graphics::arrows(x0, y0, x1, y1, col = e_col, lwd = e_lwd,
                         length = arrow_length, angle = 20)
      } else {
        graphics::segments(x0, y0, x1, y1, col = e_col, lwd = e_lwd)
      }
    }
  }

  if (n > 0L) {
    hal_idx <- which(v_vis & v_hal)
    if (length(hal_idx) > 0L) {
      graphics::points(pos[hal_idx, 1L], pos[hal_idx, 2L],
                       pch = 19L, col = v_hcol[hal_idx], cex = v_hsiz[hal_idx])
    }
    vis_idx <- which(v_vis)
    if (length(vis_idx) > 0L) {
      filled <- v_shp[vis_idx] %in% 21:25
      if (any(filled)) {
        fi <- vis_idx[filled]
        graphics::points(pos[fi, 1L], pos[fi, 2L],
                         pch = v_shp[fi], bg = v_col[fi], col = v_brd[fi],
                         cex = v_siz[fi], ...)
      }
      if (any(!filled)) {
        ui <- vis_idx[!filled]
        graphics::points(pos[ui, 1L], pos[ui, 2L],
                         pch = v_shp[ui], col = v_col[ui],
                         cex = v_siz[ui], ...)
      }
      if (!is.null(vertex_label)) {
        v_lab <- rep_len(as.character(vertex_label), n)
        graphics::text(pos[vis_idx, 1L], pos[vis_idx, 2L],
                       labels = v_lab[vis_idx], cex = vertex_label_cex)
      }
    }
  }

  if (!is.null(main) && !add) {
    graphics::title(main = main)
  }
  invisible(pos)
}

build_sirs_html <- function(pos, edges, states, newly_infected, fear, counts,
                            hide_recovered, fear_state, width, height, fps, main) {
  n <- nrow(pos)
  n_frames <- ncol(states)
  frames_s <- character(n_frames)
  frames_n <- character(n_frames)
  frames_f <- character(n_frames)
  for (t in seq_len(n_frames)) {
    frames_s[t] <- json_int_vec(states[, t])
    frames_n[t] <- json_int_vec(which(newly_infected[, t]) - 1L)
    frames_f[t] <- json_int_vec(which(fear[, t]) - 1L)
  }
  data_json <- paste0(
    "{",
    "\"pos\":", json_pos_matrix(pos), ",",
    "\"edges\":", json_edge_matrix(edges), ",",
    "\"states\":[", paste(frames_s, collapse = ","), "],",
    "\"newly\":[", paste(frames_n, collapse = ","), "],",
    "\"fear\":[", paste(frames_f, collapse = ","), "],",
    "\"countsS\":", json_int_vec(counts[, 1L]), ",",
    "\"countsI\":", json_int_vec(counts[, 2L]), ",",
    "\"countsR\":", json_int_vec(counts[, 3L]), ",",
    "\"hideRecovered\":", if (isTRUE(hide_recovered)) "true" else "false", ",",
    "\"fearMode\":", if (isTRUE(fear_state)) "true" else "false", ",",
    "\"fps\":", sprintf("%.3g", fps), ",",
    "\"n\":", as.integer(n),
    "}"
  )

  title_esc <- html_escape(main)
  paste(c(
    "<!DOCTYPE html>",
    "<html lang=\"en\">",
    "<head>",
    "<meta charset=\"utf-8\">",
    paste0("<title>", title_esc, "</title>"),
    "<style>",
    "body{margin:0;padding:18px;font-family:system-ui,-apple-system,sans-serif;background:#f8fafc;color:#0f172a;}",
    ".card{max-width:960px;margin:0 auto;background:#fff;border:1px solid #e2e8f0;border-radius:12px;padding:16px 20px;box-shadow:0 2px 8px rgba(15,23,42,.06);}",
    "h1{margin:0 0 10px;font-size:18px;font-weight:600;}",
    ".controls{display:flex;flex-wrap:wrap;gap:10px;align-items:center;margin:10px 0 12px;font-size:13px;}",
    "button{padding:6px 12px;border:1px solid #cbd5e1;border-radius:6px;background:#fff;cursor:pointer;font-size:13px;}",
    "button:hover{background:#f1f5f9;}",
    ".legend{display:flex;gap:14px;align-items:center;font-size:13px;margin-bottom:8px;flex-wrap:wrap;}",
    ".dot{display:inline-block;width:10px;height:10px;border-radius:50%;margin-right:5px;vertical-align:middle;}",
    "canvas{display:block;width:100%;border:1px solid #e2e8f0;border-radius:8px;background:#ffffff;}",
    ".spark{margin-top:10px;}",
    "</style>",
    "</head>",
    "<body>",
    "<div class=\"card\">",
    paste0("<h1>", title_esc, "</h1>"),
    "<div class=\"legend\">",
    "<span><span class=\"dot\" style=\"background:#3b82f6;\"></span>Susceptible (S): <strong id=\"cS\">0</strong></span>",
    "<span><span class=\"dot\" style=\"background:#111827;\"></span>Infected (I): <strong id=\"cI\">0</strong></span>",
    "<span><span class=\"dot\" style=\"background:#cbd5e1;border:1px solid #94a3b8;\"></span>Recovered (R): <strong id=\"cR\">0</strong></span>",
    "<span><span class=\"dot\" style=\"background:rgba(220,38,38,0.65);\"></span>Newly infected halo</span>",
    if (isTRUE(fear_state)) "<span><span class=\"dot\" style=\"background:#f59e0b;\"></span>Exposed/Fear (S near I)</span>" else "",
    "</div>",
    paste0("<canvas id=\"cv\" width=\"", as.integer(width), "\" height=\"", as.integer(height), "\"></canvas>"),
    paste0("<canvas id=\"sp\" class=\"spark\" width=\"", as.integer(width), "\" height=\"90\"></canvas>"),
    "<div class=\"controls\">",
    "<button id=\"playBtn\">Pause</button>",
    "<button id=\"resetBtn\">Reset</button>",
    "<label>Step <input id=\"stepRange\" type=\"range\" min=\"0\" max=\"0\" value=\"0\" style=\"width:220px;\"> <span id=\"stepLbl\">0</span></label>",
    "<label>FPS <input id=\"fpsRange\" type=\"range\" min=\"1\" max=\"30\" value=\"12\" style=\"width:90px;\"> <span id=\"fpsLbl\">12</span></label>",
    "<label><input id=\"hideChk\" type=\"checkbox\"> Hide recovered</label>",
    "</div>",
    "</div>",
    "<script>",
    paste0("const D = ", data_json, ";"),
    "const cv = document.getElementById('cv'), ctx = cv.getContext('2d');",
    "const sp = document.getElementById('sp'), sctx = sp.getContext('2d');",
    "const stepRange = document.getElementById('stepRange'), stepLbl = document.getElementById('stepLbl');",
    "const fpsRange = document.getElementById('fpsRange'), fpsLbl = document.getElementById('fpsLbl');",
    "const hideChk = document.getElementById('hideChk'), playBtn = document.getElementById('playBtn');",
    "const nFrames = D.states.length;",
    "stepRange.max = Math.max(0, nFrames - 1);",
    "fpsRange.value = Math.max(1, Math.min(30, Math.round(D.fps)));",
    "fpsLbl.textContent = fpsRange.value;",
    "hideChk.checked = D.hideRecovered;",
    "let frame = 0, playing = true, lastTime = 0;",
    "let minX=Infinity, maxX=-Infinity, minY=Infinity, maxY=-Infinity;",
    "for (const p of D.pos){if(p[0]<minX)minX=p[0];if(p[0]>maxX)maxX=p[0];if(p[1]<minY)minY=p[1];if(p[1]>maxY)maxY=p[1];}",
    "if(!isFinite(minX)){minX=-1;maxX=1;minY=-1;maxY=1;}",
    "const spanX = Math.max(1e-6, maxX - minX), spanY = Math.max(1e-6, maxY - minY);",
    "const pad = 34, scale = Math.min((cv.width - 2*pad)/spanX, (cv.height - 2*pad)/spanY);",
    "const ox = (cv.width - spanX*scale)/2 - minX*scale, oy = (cv.height - spanY*scale)/2 - minY*scale;",
    "const pts = D.pos.map(p => [p[0]*scale + ox, cv.height - (p[1]*scale + oy)]);",
    "const rNode = Math.max(3, Math.min(7, 55 / Math.sqrt(Math.max(1, D.n))));",
    "function drawSpark(t){",
    "  sctx.clearRect(0,0,sp.width,sp.height);",
    "  const m=10, w=sp.width-2*m, h=sp.height-2*m, N=Math.max(1,D.n);",
    "  const series=[[D.countsS,'#3b82f6'],[D.countsI,'#111827'],[D.countsR,'#94a3b8']];",
    "  for(const [arr,col] of series){",
    "    sctx.beginPath(); sctx.strokeStyle=col; sctx.lineWidth=2;",
    "    for(let i=0;i<nFrames;i++){",
    "      const x = m + (nFrames>1 ? i/(nFrames-1)*w : 0);",
    "      const y = m + h*(1 - arr[i]/N);",
    "      if(i===0) sctx.moveTo(x,y); else sctx.lineTo(x,y);",
    "    }",
    "    sctx.stroke();",
    "  }",
    "  const cx = m + (nFrames>1 ? t/(nFrames-1)*w : 0);",
    "  sctx.strokeStyle='#dc2626'; sctx.lineWidth=1.5; sctx.beginPath(); sctx.moveTo(cx,m); sctx.lineTo(cx,m+h); sctx.stroke();",
    "}",
    "function render(t){",
    "  ctx.clearRect(0,0,cv.width,cv.height);",
    "  const st = D.states[t];",
    "  const newlySet = new Set(D.newly[t]);",
    "  const fearSet = new Set(D.fear[t]);",
    "  const hideR = hideChk.checked;",
    "  ctx.strokeStyle = 'rgba(100,116,139,0.28)'; ctx.lineWidth = 1;",
    "  ctx.beginPath();",
    "  for(const e of D.edges){",
    "    const u=e[0], v=e[1];",
    "    if(hideR && (st[u]===2 || st[v]===2)) continue;",
    "    ctx.moveTo(pts[u][0], pts[u][1]); ctx.lineTo(pts[v][0], pts[v][1]);",
    "  }",
    "  ctx.stroke();",
    "  for(const u of newlySet){",
    "    if(hideR && st[u]===2) continue;",
    "    ctx.fillStyle='rgba(220,38,38,0.35)';",
    "    ctx.beginPath(); ctx.arc(pts[u][0], pts[u][1], rNode*2.2, 0, Math.PI*2); ctx.fill();",
    "  }",
    "  for(let u=0;u<D.n;u++){",
    "    if(hideR && st[u]===2) continue;",
    "    let fill = '#3b82f6';",
    "    if(st[u]===1) fill = '#111827';",
    "    else if(st[u]===2) fill = '#cbd5e1';",
    "    else if(D.fearMode && fearSet.has(u)) fill = '#f59e0b';",
    "    ctx.fillStyle = fill;",
    "    ctx.beginPath(); ctx.arc(pts[u][0], pts[u][1], rNode, 0, Math.PI*2); ctx.fill();",
    "    ctx.strokeStyle = '#1e293b'; ctx.lineWidth = 0.8; ctx.stroke();",
    "  }",
    "  document.getElementById('cS').textContent = D.countsS[t];",
    "  document.getElementById('cI').textContent = D.countsI[t];",
    "  document.getElementById('cR').textContent = D.countsR[t];",
    "  stepRange.value = t; stepLbl.textContent = t + ' / ' + (nFrames - 1);",
    "  drawSpark(t);",
    "}",
    "function tick(ts){",
    "  if(playing && nFrames > 1){",
    "    const interval = 1000 / Number(fpsRange.value);",
    "    if(ts - lastTime >= interval){ lastTime = ts; frame = (frame + 1) % nFrames; render(frame); }",
    "  }",
    "  requestAnimationFrame(tick);",
    "}",
    "playBtn.onclick = () => { playing = !playing; playBtn.textContent = playing ? 'Pause' : 'Play'; };",
    "document.getElementById('resetBtn').onclick = () => { frame = 0; render(0); };",
    "stepRange.oninput = () => { playing = false; playBtn.textContent = 'Play'; frame = Number(stepRange.value); render(frame); };",
    "fpsRange.oninput = () => { fpsLbl.textContent = fpsRange.value; };",
    "hideChk.onchange = () => { render(frame); };",
    "render(0); requestAnimationFrame(tick);",
    "</script>",
    "</body>",
    "</html>"
  ), collapse = "\n")
}

animate_sirs <- function(graph, pos = NULL, n_steps = 120L,
                         x = 0.001, r = 0.1, s = 0.01,
                         init_state = NULL, directed = FALSE,
                         hide_recovered = TRUE, fear_state = FALSE,
                         mode = c("html", "frames", "device", "none"),
                         output = NULL, width = 760L, height = 560L,
                         fps = 12, delay = 0.05, open = FALSE,
                         main = "SIRS Epidemic Animation") {
  graph_pointer(graph)
  mode <- match.arg(mode)
  n <- graph_num_vertices(graph)
  n_steps <- validate_nonneg_int(n_steps, "n_steps", 0L)
  x <- validate_prob_scalar(x, "x")
  r <- validate_prob_scalar(r, "r")
  s <- validate_prob_scalar(s, "s")
  width <- validate_nonneg_int(width, "width", 100L)
  height <- validate_nonneg_int(height, "height", 100L)
  fps <- validate_pos_scalar(fps, "fps")
  delay <- validate_pos_scalar(delay, "delay", allow_zero = TRUE)
  if (!is.logical(directed) || length(directed) != 1L || is.na(directed)) {
    stop("directed must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(hide_recovered) || length(hide_recovered) != 1L || is.na(hide_recovered)) {
    stop("hide_recovered must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(fear_state) || length(fear_state) != 1L || is.na(fear_state)) {
    stop("fear_state must be TRUE or FALSE", call. = FALSE)
  }
  if (is.null(pos)) {
    pos <- sfdp_layout(graph)
  } else {
    pos <- validate_pos_matrix(pos, n, "pos")
  }
  if (is.null(init_state)) {
    init_state <- integer(n)
  } else {
    if (!is.numeric(init_state) || length(init_state) != n || anyNA(init_state) ||
        any(!(init_state %in% c(0, 1, 2)))) {
      stop("init_state must be an integer vector of length N with values in {0, 1, 2}", call. = FALSE)
    }
    init_state <- as.integer(init_state)
  }

  sim <- cpp_sirs_animation_sim(
    pointer = graph$pointer,
    n_steps = n_steps,
    x_spont = x,
    r_rec = r,
    s_resusc = s,
    init_state = init_state,
    directed = directed
  )

  edges <- graph_edges(graph)
  files <- character(0)
  out_path <- output
  lims <- compute_plot_limits(pos, pad = 0.08)

  draw_sirs_frame <- function(t) {
    st <- sim$states[, t]
    v_col <- ifelse(st == 1L, "#111827", ifelse(st == 2L, "#CBD5E1", "#3B82F6"))
    if (isTRUE(fear_state)) {
      v_col[sim$fear[, t] & st == 0L] <- "#F59E0B"
    }
    v_vis <- if (isTRUE(hide_recovered)) (st != 2L) else rep(TRUE, n)
    plot_graph(
      graph, pos = pos,
      vertex_color = v_col,
      vertex_halo = sim$newly_infected[, t],
      vertex_halo_color = grDevices::adjustcolor("#DC2626", 0.40),
      vertex_halo_size = 2.3,
      vertex_visible = v_vis,
      xlim = lims$xlim, ylim = lims$ylim,
      main = sprintf("%s (step %d: S=%d, I=%d, R=%d)",
                     main, t - 1L, sim$counts[t, 1L], sim$counts[t, 2L], sim$counts[t, 3L])
    )
  }

  if (mode == "html") {
    if (is.null(out_path)) {
      out_path <- tempfile("rgraphtool_sirs_", fileext = ".html")
    }
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    html_txt <- build_sirs_html(
      pos = pos, edges = edges,
      states = sim$states, newly_infected = sim$newly_infected,
      fear = sim$fear, counts = sim$counts,
      hide_recovered = hide_recovered, fear_state = fear_state,
      width = width, height = height, fps = fps, main = main
    )
    writeLines(html_txt, con = out_path, useBytes = TRUE)
    files <- out_path
    if (isTRUE(open) && interactive()) {
      utils::browseURL(out_path)
    }
  } else if (mode == "frames") {
    if (is.null(out_path)) {
      out_path <- tempfile("rgraphtool_sirs_frames_")
    }
    dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
    n_frames <- ncol(sim$states)
    files <- file.path(out_path, sprintf("sirs_%06d.png", seq_len(n_frames) - 1L))
    for (t in seq_len(n_frames)) {
      grDevices::png(filename = files[t], width = width, height = height)
      tryCatch(draw_sirs_frame(t), finally = grDevices::dev.off())
    }
  } else if (mode == "device") {
    n_frames <- ncol(sim$states)
    for (t in seq_len(n_frames)) {
      grDevices::dev.hold()
      draw_sirs_frame(t)
      grDevices::dev.flush()
      if (delay > 0 && t < n_frames) Sys.sleep(delay)
    }
  }

  structure(
    list(
      pos = pos,
      states = sim$states,
      newly_infected = sim$newly_infected,
      fear = sim$fear,
      counts = sim$counts,
      params = list(n_steps = n_steps, x = x, r = r, s = s,
                    directed = directed, hide_recovered = hide_recovered,
                    fear_state = fear_state),
      mode = mode,
      output = out_path,
      files = files
    ),
    class = "rgraphtool_sirs_animation"
  )
}

print.rgraphtool_sirs_animation <- function(x, ...) {
  n <- nrow(x$states)
  nf <- ncol(x$states)
  peak_i <- max(x$counts[, "I"])
  total_new <- sum(x$newly_infected)
  cat("rgraphtool SIRS epidemic animation\n")
  cat(sprintf("  vertices: %d | frames: %d (steps = %d)\n", n, nf, x$params$n_steps))
  cat(sprintf("  rates: x = %g, r = %g, s = %g | hide_recovered = %s | fear_state = %s\n",
              x$params$x, x$params$r, x$params$s,
              x$params$hide_recovered, x$params$fear_state))
  cat(sprintf("  peak infected: %d | total new infections: %d\n", peak_i, total_new))
  cat(sprintf("  mode: %s%s\n", x$mode,
              if (!is.null(x$output)) paste0(" -> ", x$output) else ""))
  invisible(x)
}

build_dynamic_layout_html <- function(positions, edges_list, mean_edge_length,
                                      rewired, width, height, fps, main) {
  n_frames <- length(positions)
  n <- if (n_frames > 0L) nrow(positions[[1L]]) else 0L
  pos_json <- paste0("[", paste(vapply(positions, json_pos_matrix, character(1L)), collapse = ","), "]")
  edg_json <- paste0("[", paste(vapply(edges_list, json_edge_matrix, character(1L)), collapse = ","), "]")
  data_json <- paste0(
    "{",
    "\"positions\":", pos_json, ",",
    "\"edges\":", edg_json, ",",
    "\"meanLen\":", json_num_vec(mean_edge_length, 5L), ",",
    "\"rewired\":", json_int_vec(rewired), ",",
    "\"fps\":", sprintf("%.3g", fps), ",",
    "\"n\":", as.integer(n),
    "}"
  )
  title_esc <- html_escape(main)
  paste(c(
    "<!DOCTYPE html>",
    "<html lang=\"en\">",
    "<head>",
    "<meta charset=\"utf-8\">",
    paste0("<title>", title_esc, "</title>"),
    "<style>",
    "body{margin:0;padding:18px;font-family:system-ui,-apple-system,sans-serif;background:#f8fafc;color:#0f172a;}",
    ".card{max-width:960px;margin:0 auto;background:#fff;border:1px solid #e2e8f0;border-radius:12px;padding:16px 20px;box-shadow:0 2px 8px rgba(15,23,42,.06);}",
    "h1{margin:0 0 8px;font-size:18px;font-weight:600;}",
    ".stats{display:flex;gap:18px;font-size:13px;margin-bottom:10px;color:#334155;}",
    ".controls{display:flex;flex-wrap:wrap;gap:10px;align-items:center;margin-top:12px;font-size:13px;}",
    "button{padding:6px 12px;border:1px solid #cbd5e1;border-radius:6px;background:#fff;cursor:pointer;font-size:13px;}",
    "button:hover{background:#f1f5f9;}",
    "canvas{display:block;width:100%;border:1px solid #e2e8f0;border-radius:8px;background:#ffffff;}",
    ".spark{margin-top:10px;}",
    "</style>",
    "</head>",
    "<body>",
    "<div class=\"card\">",
    paste0("<h1>", title_esc, "</h1>"),
    "<div class=\"stats\">",
    "<span>Mean edge length: <strong id=\"lenLbl\">0</strong></span>",
    "<span>Cumulative rewirings: <strong id=\"rewLbl\">0</strong></span>",
    "</div>",
    paste0("<canvas id=\"cv\" width=\"", as.integer(width), "\" height=\"", as.integer(height), "\"></canvas>"),
    paste0("<canvas id=\"sp\" class=\"spark\" width=\"", as.integer(width), "\" height=\"85\"></canvas>"),
    "<div class=\"controls\">",
    "<button id=\"playBtn\">Pause</button>",
    "<button id=\"resetBtn\">Reset</button>",
    "<label>Step <input id=\"stepRange\" type=\"range\" min=\"0\" max=\"0\" value=\"0\" style=\"width:220px;\"> <span id=\"stepLbl\">0</span></label>",
    "<label>FPS <input id=\"fpsRange\" type=\"range\" min=\"1\" max=\"30\" value=\"15\" style=\"width:90px;\"> <span id=\"fpsLbl\">15</span></label>",
    "<label><input id=\"autofitChk\" type=\"checkbox\" checked> Auto-fit frame</label>",
    "</div>",
    "</div>",
    "<script>",
    paste0("const D = ", data_json, ";"),
    "const cv = document.getElementById('cv'), ctx = cv.getContext('2d');",
    "const sp = document.getElementById('sp'), sctx = sp.getContext('2d');",
    "const stepRange = document.getElementById('stepRange'), stepLbl = document.getElementById('stepLbl');",
    "const fpsRange = document.getElementById('fpsRange'), fpsLbl = document.getElementById('fpsLbl');",
    "const autofitChk = document.getElementById('autofitChk'), playBtn = document.getElementById('playBtn');",
    "const nFrames = D.positions.length;",
    "stepRange.max = Math.max(0, nFrames - 1);",
    "fpsRange.value = Math.max(1, Math.min(30, Math.round(D.fps)));",
    "fpsLbl.textContent = fpsRange.value;",
    "const cumRew = []; let acc = 0; for (const r of D.rewired){ acc += r; cumRew.push(acc); }",
    "let gMinX=Infinity,gMaxX=-Infinity,gMinY=Infinity,gMaxY=-Infinity;",
    "for(const fr of D.positions) for(const p of fr){",
    "  if(p[0]<gMinX)gMinX=p[0]; if(p[0]>gMaxX)gMaxX=p[0];",
    "  if(p[1]<gMinY)gMinY=p[1]; if(p[1]>gMaxY)gMaxY=p[1];",
    "}",
    "if(!isFinite(gMinX)){gMinX=-1;gMaxX=1;gMinY=-1;gMaxY=1;}",
    "const maxLen = Math.max(1e-6, ...D.meanLen), minLen = Math.min(...D.meanLen);",
    "const rNode = Math.max(3, Math.min(6.5, 50 / Math.sqrt(Math.max(1, D.n))));",
    "let frame = 0, playing = true, lastTime = 0;",
    "function drawSpark(t){",
    "  sctx.clearRect(0,0,sp.width,sp.height);",
    "  const m=10, w=sp.width-2*m, h=sp.height-2*m, span=Math.max(1e-6, maxLen - Math.min(0, minLen));",
    "  sctx.beginPath(); sctx.strokeStyle='#2563eb'; sctx.lineWidth=2;",
    "  for(let i=0;i<nFrames;i++){",
    "    const x = m + (nFrames>1 ? i/(nFrames-1)*w : 0);",
    "    const y = m + h*(1 - (D.meanLen[i] - Math.min(0, minLen))/span);",
    "    if(i===0) sctx.moveTo(x,y); else sctx.lineTo(x,y);",
    "  }",
    "  sctx.stroke();",
    "  const cx = m + (nFrames>1 ? t/(nFrames-1)*w : 0);",
    "  sctx.strokeStyle='#dc2626'; sctx.lineWidth=1.5; sctx.beginPath(); sctx.moveTo(cx,m); sctx.lineTo(cx,m+h); sctx.stroke();",
    "}",
    "function render(t){",
    "  ctx.clearRect(0,0,cv.width,cv.height);",
    "  const pos = D.positions[t], edg = D.edges[t];",
    "  let minX=gMinX, maxX=gMaxX, minY=gMinY, maxY=gMaxY;",
    "  if(autofitChk.checked && pos.length > 0){",
    "    minX=Infinity; maxX=-Infinity; minY=Infinity; maxY=-Infinity;",
    "    for(const p of pos){if(p[0]<minX)minX=p[0];if(p[0]>maxX)maxX=p[0];if(p[1]<minY)minY=p[1];if(p[1]>maxY)maxY=p[1];}",
    "  }",
    "  const spanX = Math.max(1e-6, maxX - minX), spanY = Math.max(1e-6, maxY - minY);",
    "  const pad = 32, scale = Math.min((cv.width - 2*pad)/spanX, (cv.height - 2*pad)/spanY);",
    "  const ox = (cv.width - spanX*scale)/2 - minX*scale, oy = (cv.height - spanY*scale)/2 - minY*scale;",
    "  const pts = pos.map(p => [p[0]*scale + ox, cv.height - (p[1]*scale + oy)]);",
    "  ctx.strokeStyle='rgba(71,85,105,0.35)'; ctx.lineWidth=1.1; ctx.beginPath();",
    "  for(const e of edg){ const a=pts[e[0]], b=pts[e[1]]; ctx.moveTo(a[0],a[1]); ctx.lineTo(b[0],b[1]); }",
    "  ctx.stroke();",
    "  for(let u=0;u<pts.length;u++){",
    "    ctx.fillStyle='#2563eb'; ctx.beginPath(); ctx.arc(pts[u][0],pts[u][1],rNode,0,Math.PI*2); ctx.fill();",
    "    ctx.strokeStyle='#1e293b'; ctx.lineWidth=0.7; ctx.stroke();",
    "  }",
    "  document.getElementById('lenLbl').textContent = D.meanLen[t].toFixed(3);",
    "  document.getElementById('rewLbl').textContent = cumRew[t];",
    "  stepRange.value = t; stepLbl.textContent = t + ' / ' + (nFrames - 1);",
    "  drawSpark(t);",
    "}",
    "function tick(ts){",
    "  if(playing && nFrames > 1){",
    "    const interval = 1000 / Number(fpsRange.value);",
    "    if(ts - lastTime >= interval){ lastTime = ts; frame = (frame + 1) % nFrames; render(frame); }",
    "  }",
    "  requestAnimationFrame(tick);",
    "}",
    "playBtn.onclick = () => { playing = !playing; playBtn.textContent = playing ? 'Pause' : 'Play'; };",
    "document.getElementById('resetBtn').onclick = () => { frame = 0; render(0); };",
    "stepRange.oninput = () => { playing = false; playBtn.textContent = 'Play'; frame = Number(stepRange.value); render(frame); };",
    "fpsRange.oninput = () => { fpsLbl.textContent = fpsRange.value; };",
    "autofitChk.onchange = () => { render(frame); };",
    "render(0); requestAnimationFrame(tick);",
    "</script>",
    "</body>",
    "</html>"
  ), collapse = "\n")
}

animate_dynamic_layout <- function(graph, pos = NULL, n_steps = 80L,
                                   K = 0.5, C = 0.2, p = 2.0,
                                   step = 0.005, n_rewire = 10L,
                                   mode = c("html", "frames", "device", "none"),
                                   output = NULL, width = 760L, height = 560L,
                                   fps = 15, delay = 0.05, open = FALSE,
                                   main = "Dynamic Layout & Spatial Segregation") {
  graph_pointer(graph)
  mode <- match.arg(mode)
  n <- graph_num_vertices(graph)
  n_steps <- validate_nonneg_int(n_steps, "n_steps", 0L)
  K <- validate_pos_scalar(K, "K")
  C <- validate_pos_scalar(C, "C")
  p <- validate_pos_scalar(p, "p")
  step <- validate_pos_scalar(step, "step")
  n_rewire <- validate_nonneg_int(n_rewire, "n_rewire", 0L)
  width <- validate_nonneg_int(width, "width", 100L)
  height <- validate_nonneg_int(height, "height", 100L)
  fps <- validate_pos_scalar(fps, "fps")
  delay <- validate_pos_scalar(delay, "delay", allow_zero = TRUE)

  if (is.null(pos)) {
    pos <- sfdp_layout(graph, K = K, C = C, p = p)
  } else {
    pos <- validate_pos_matrix(pos, n, "pos")
  }

  sim <- cpp_dynamic_layout_sim(
    pointer = graph$pointer,
    init_pos = pos,
    n_steps = n_steps,
    K = K,
    C = C,
    p = p,
    step = step,
    n_rewire = n_rewire
  )

  draw_dyn_frame <- function(t) {
    pm <- sim$positions[[t]]
    em <- sim$edges[[t]]
    lims <- compute_plot_limits(pm, pad = 0.08)
    op <- graphics::par(bg = "white", mar = c(1, 1, 2.5, 1))
    on.exit(graphics::par(op), add = TRUE)
    graphics::plot.new()
    graphics::plot.window(xlim = lims$xlim, ylim = lims$ylim, asp = 1)
    if (nrow(em) > 0L) {
      graphics::segments(pm[em[, 1L], 1L], pm[em[, 1L], 2L],
                         pm[em[, 2L], 1L], pm[em[, 2L], 2L],
                         col = grDevices::adjustcolor("#475569", 0.35), lwd = 1.0)
    }
    if (nrow(pm) > 0L) {
      graphics::points(pm[, 1L], pm[, 2L], pch = 21L,
                       bg = "#2563EB", col = "#1E293B", cex = 1.1)
    }
    graphics::title(main = sprintf("%s (step %d: mean edge len = %.3f)",
                                   main, t - 1L, sim$mean_edge_length[t]))
  }

  files <- character(0)
  out_path <- output
  n_frames <- length(sim$positions)

  if (mode == "html") {
    if (is.null(out_path)) {
      out_path <- tempfile("rgraphtool_dynlayout_", fileext = ".html")
    }
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    html_txt <- build_dynamic_layout_html(
      positions = sim$positions,
      edges_list = sim$edges,
      mean_edge_length = sim$mean_edge_length,
      rewired = sim$rewired,
      width = width, height = height, fps = fps, main = main
    )
    writeLines(html_txt, con = out_path, useBytes = TRUE)
    files <- out_path
    if (isTRUE(open) && interactive()) {
      utils::browseURL(out_path)
    }
  } else if (mode == "frames") {
    if (is.null(out_path)) {
      out_path <- tempfile("rgraphtool_dynlayout_frames_")
    }
    dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
    files <- file.path(out_path, sprintf("dance_%06d.png", seq_len(n_frames) - 1L))
    for (t in seq_len(n_frames)) {
      grDevices::png(filename = files[t], width = width, height = height)
      tryCatch(draw_dyn_frame(t), finally = grDevices::dev.off())
    }
  } else if (mode == "device") {
    for (t in seq_len(n_frames)) {
      grDevices::dev.hold()
      draw_dyn_frame(t)
      grDevices::dev.flush()
      if (delay > 0 && t < n_frames) Sys.sleep(delay)
    }
  }

  structure(
    list(
      positions = sim$positions,
      edges = sim$edges,
      mean_edge_length = sim$mean_edge_length,
      rewired = sim$rewired,
      params = list(n_steps = n_steps, K = K, C = C, p = p,
                    step = step, n_rewire = n_rewire),
      mode = mode,
      output = out_path,
      files = files
    ),
    class = "rgraphtool_dynamic_layout"
  )
}

print.rgraphtool_dynamic_layout <- function(x, ...) {
  nf <- length(x$positions)
  n <- if (nf > 0L) nrow(x$positions[[1L]]) else 0L
  m <- if (nf > 0L) nrow(x$edges[[1L]]) else 0L
  cat("rgraphtool dynamic layout animation\n")
  cat(sprintf("  vertices: %d | edges: %d | frames: %d\n", n, m, nf))
  if (nf > 0L) {
    cat(sprintf("  mean edge length: %.4f -> %.4f | total rewirings: %d\n",
                x$mean_edge_length[1L], x$mean_edge_length[nf], sum(x$rewired)))
  }
  cat(sprintf("  mode: %s%s\n", x$mode,
              if (!is.null(x$output)) paste0(" -> ", x$output) else ""))
  invisible(x)
}

build_interactive_bfs_html <- function(pos, edges, source, max_dist, directed,
                                       highlight_color, base_color,
                                       width, height, main) {
  n <- nrow(pos)
  md_init <- if (is.infinite(max_dist)) max(1L, n) else as.integer(max_dist)
  data_json <- paste0(
    "{",
    "\"pos\":", json_pos_matrix(pos), ",",
    "\"edges\":", json_edge_matrix(edges), ",",
    "\"source\":", as.integer(source - 1L), ",",
    "\"maxDist\":", as.integer(md_init), ",",
    "\"directed\":", if (isTRUE(directed)) "true" else "false", ",",
    "\"hiColor\":\"", html_escape(highlight_color), "\",",
    "\"baseColor\":\"", html_escape(base_color), "\",",
    "\"n\":", as.integer(n),
    "}"
  )
  title_esc <- html_escape(main)
  paste(c(
    "<!DOCTYPE html>",
    "<html lang=\"en\">",
    "<head>",
    "<meta charset=\"utf-8\">",
    paste0("<title>", title_esc, "</title>"),
    "<style>",
    "body{margin:0;padding:18px;font-family:system-ui,-apple-system,sans-serif;background:#f8fafc;color:#0f172a;}",
    ".card{max-width:960px;margin:0 auto;background:#fff;border:1px solid #e2e8f0;border-radius:12px;padding:16px 20px;box-shadow:0 2px 8px rgba(15,23,42,.06);}",
    "h1{margin:0 0 8px;font-size:18px;font-weight:600;}",
    ".info{display:flex;flex-wrap:wrap;gap:16px;font-size:13px;margin-bottom:10px;color:#334155;}",
    ".controls{display:flex;flex-wrap:wrap;gap:14px;align-items:center;margin-top:10px;font-size:13px;}",
    "canvas{display:block;width:100%;border:1px solid #e2e8f0;border-radius:8px;background:#ffffff;cursor:crosshair;}",
    "</style>",
    "</head>",
    "<body>",
    "<div class=\"card\">",
    paste0("<h1>", title_esc, "</h1>"),
    "<div class=\"info\">",
    "<span>Hover or click any vertex to highlight its BFS tree.</span>",
    "<span>Source vertex: <strong id=\"srcLbl\">1</strong></span>",
    "<span>Reached: <strong id=\"reachLbl\">0</strong></span>",
    "<span>Pinned: <strong id=\"pinLbl\">No</strong> (click to toggle)</span>",
    "</div>",
    paste0("<canvas id=\"cv\" width=\"", as.integer(width), "\" height=\"", as.integer(height), "\"></canvas>"),
    "<div class=\"controls\">",
    "<label>Max BFS depth: <input id=\"distRange\" type=\"range\" min=\"1\" max=\"12\" value=\"3\" style=\"width:180px;\"> <span id=\"distLbl\">3</span></label>",
    "</div>",
    "</div>",
    "<script>",
    paste0("const D = ", data_json, ";"),
    "const cv = document.getElementById('cv'), ctx = cv.getContext('2d');",
    "const distRange = document.getElementById('distRange'), distLbl = document.getElementById('distLbl');",
    "distRange.max = Math.max(10, Math.min(25, D.maxDist + 4));",
    "distRange.value = Math.min(Number(distRange.max), Math.max(1, D.maxDist));",
    "distLbl.textContent = distRange.value;",
    "const adj = Array.from({length: D.n}, () => []);",
    "for (let i = 0; i < D.edges.length; i++){",
    "  const u = D.edges[i][0], v = D.edges[i][1];",
    "  if (u === v) continue;",
    "  adj[u].push([v, i]);",
    "  if (!D.directed) adj[v].push([u, i]);",
    "}",
    "let minX=Infinity, maxX=-Infinity, minY=Infinity, maxY=-Infinity;",
    "for (const p of D.pos){if(p[0]<minX)minX=p[0];if(p[0]>maxX)maxX=p[0];if(p[1]<minY)minY=p[1];if(p[1]>maxY)maxY=p[1];}",
    "if(!isFinite(minX)){minX=-1;maxX=1;minY=-1;maxY=1;}",
    "const spanX = Math.max(1e-6, maxX - minX), spanY = Math.max(1e-6, maxY - minY);",
    "const pad = 34, scale = Math.min((cv.width - 2*pad)/spanX, (cv.height - 2*pad)/spanY);",
    "const ox = (cv.width - spanX*scale)/2 - minX*scale, oy = (cv.height - spanY*scale)/2 - minY*scale;",
    "const pts = D.pos.map(p => [p[0]*scale + ox, cv.height - (p[1]*scale + oy)]);",
    "const rNode = Math.max(3.5, Math.min(7.5, 55 / Math.sqrt(Math.max(1, D.n))));",
    "let currentSrc = D.source, pinned = false;",
    "function runBfs(src, maxD){",
    "  const dist = new Int32Array(D.n).fill(-1);",
    "  const treeEdge = new Uint8Array(D.edges.length);",
    "  const q = [src]; dist[src] = 0; let head = 0;",
    "  while(head < q.length){",
    "    const u = q[head++];",
    "    if(dist[u] >= maxD) continue;",
    "    for(const [v, ei] of adj[u]){",
    "      if(dist[v] === -1){ dist[v] = dist[u] + 1; treeEdge[ei] = 1; q.push(v); }",
    "    }",
    "  }",
    "  return {dist, treeEdge, reached: q.length};",
    "}",
    "function render(src){",
    "  const maxD = Number(distRange.value);",
    "  const res = runBfs(src, maxD);",
    "  ctx.clearRect(0,0,cv.width,cv.height);",
    "  ctx.strokeStyle = 'rgba(153,153,153,0.45)'; ctx.lineWidth = 1; ctx.beginPath();",
    "  for(let i=0;i<D.edges.length;i++){",
    "    if(res.treeEdge[i]) continue;",
    "    const a=pts[D.edges[i][0]], b=pts[D.edges[i][1]];",
    "    ctx.moveTo(a[0],a[1]); ctx.lineTo(b[0],b[1]);",
    "  }",
    "  ctx.stroke();",
    "  ctx.strokeStyle = D.hiColor; ctx.lineWidth = 2.6; ctx.beginPath();",
    "  for(let i=0;i<D.edges.length;i++){",
    "    if(!res.treeEdge[i]) continue;",
    "    const a=pts[D.edges[i][0]], b=pts[D.edges[i][1]];",
    "    ctx.moveTo(a[0],a[1]); ctx.lineTo(b[0],b[1]);",
    "  }",
    "  ctx.stroke();",
    "  for(let u=0;u<D.n;u++){",
    "    const inTree = res.dist[u] >= 0;",
    "    ctx.fillStyle = inTree ? D.hiColor : D.baseColor;",
    "    const rad = (u === src) ? rNode * 1.65 : (inTree ? rNode * 1.15 : rNode);",
    "    ctx.beginPath(); ctx.arc(pts[u][0], pts[u][1], rad, 0, Math.PI*2); ctx.fill();",
    "    ctx.strokeStyle = (u === src) ? '#0f172a' : '#334155'; ctx.lineWidth = (u === src) ? 2 : 0.7; ctx.stroke();",
    "  }",
    "  document.getElementById('srcLbl').textContent = String(src + 1);",
    "  document.getElementById('reachLbl').textContent = res.reached + ' / ' + D.n;",
    "  document.getElementById('pinLbl').textContent = pinned ? 'Yes' : 'No';",
    "}",
    "function nearestVertex(evt){",
    "  const rect = cv.getBoundingClientRect();",
    "  const x = (evt.clientX - rect.left) * (cv.width / rect.width);",
    "  const y = (evt.clientY - rect.top) * (cv.height / rect.height);",
    "  let best = -1, bestD2 = 400;",
    "  for(let u=0;u<D.n;u++){",
    "    const dx = pts[u][0]-x, dy = pts[u][1]-y, d2 = dx*dx + dy*dy;",
    "    if(d2 < bestD2){ bestD2 = d2; best = u; }",
    "  }",
    "  return best;",
    "}",
    "cv.addEventListener('mousemove', evt => {",
    "  if(pinned) return;",
    "  const v = nearestVertex(evt);",
    "  if(v >= 0 && v !== currentSrc){ currentSrc = v; render(currentSrc); }",
    "});",
    "cv.addEventListener('click', evt => {",
    "  const v = nearestVertex(evt);",
    "  if(v >= 0){ currentSrc = v; pinned = !pinned; render(currentSrc); }",
    "  else { pinned = false; render(currentSrc); }",
    "});",
    "distRange.oninput = () => { distLbl.textContent = distRange.value; render(currentSrc); };",
    "render(currentSrc);",
    "</script>",
    "</body>",
    "</html>"
  ), collapse = "\n")
}

interactive_bfs <- function(graph, pos = NULL, source = 1L, max_dist = 3L,
                            directed = NULL,
                            mode = c("html", "locator", "snapshot"),
                            output = NULL, width = 760L, height = 560L,
                            highlight_color = "#CE5C00",
                            base_color = "#999999",
                            open = FALSE,
                            main = "Interactive BFS Tree Explorer") {
  graph_pointer(graph)
  mode <- match.arg(mode)
  n <- graph_num_vertices(graph)
  if (n == 0L) stop("graph has no vertices", call. = FALSE)
  if (!is.numeric(source) || length(source) < 1L || anyNA(source) ||
      any(source < 1) || any(source > n) || any(source != floor(source))) {
    stop(sprintf("source must contain integers between 1 and %d", n), call. = FALSE)
  }
  source <- as.integer(source)
  if (is.null(directed)) {
    directed <- isTRUE(graph$directed)
  } else if (!is.logical(directed) || length(directed) != 1L || is.na(directed)) {
    stop("directed must be TRUE, FALSE, or NULL", call. = FALSE)
  }
  width <- validate_nonneg_int(width, "width", 100L)
  height <- validate_nonneg_int(height, "height", 100L)
  if (is.null(pos)) {
    pos <- sfdp_layout(graph)
  } else {
    pos <- validate_pos_matrix(pos, n, "pos")
  }

  edges <- graph_edges(graph)
  bfs <- graph_bfs(graph, source = source[1L], max_dist = max_dist, directed = directed)
  out_path <- output

  draw_bfs_panel <- function(src_v, panel_title = NULL) {
    b <- graph_bfs(graph, source = src_v, max_dist = max_dist, directed = directed)
    in_tree <- is.finite(b$distance)
    v_col <- ifelse(in_tree, highlight_color, base_color)
    v_siz <- ifelse(seq_len(n) == src_v, 1.8, ifelse(in_tree, 1.25, 0.95))
    e_col <- ifelse(b$tree_edges, highlight_color, grDevices::adjustcolor(base_color, 0.5))
    e_lwd <- ifelse(b$tree_edges, 2.2, 0.9)
    ttl <- if (is.null(panel_title)) {
      sprintf("%s (source = %d, max_dist = %s, reached = %d)",
              main, src_v,
              if (is.infinite(max_dist)) "Inf" else as.character(max_dist),
              sum(in_tree))
    } else {
      panel_title
    }
    plot_graph(
      graph, pos = pos,
      vertex_color = v_col,
      vertex_size = v_siz,
      edge_color = e_col,
      edge_lwd = e_lwd,
      arrows = directed,
      main = ttl
    )
    b
  }

  if (mode == "html") {
    if (is.null(out_path)) {
      out_path <- tempfile("rgraphtool_bfs_", fileext = ".html")
    }
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    html_txt <- build_interactive_bfs_html(
      pos = pos, edges = edges,
      source = source[1L], max_dist = max_dist, directed = directed,
      highlight_color = highlight_color, base_color = base_color,
      width = width, height = height, main = main
    )
    writeLines(html_txt, con = out_path, useBytes = TRUE)
    if (isTRUE(open) && interactive()) {
      utils::browseURL(out_path)
    }
  } else if (mode == "snapshot") {
    render_snaps <- function() {
      k <- length(source)
      if (k > 1L) {
        nr <- floor(sqrt(k))
        nc <- ceiling(k / nr)
        op <- graphics::par(mfrow = c(nr, nc))
        on.exit(graphics::par(op), add = TRUE)
      }
      last_b <- bfs
      for (sv in source) {
        last_b <- draw_bfs_panel(sv)
      }
      last_b
    }
    if (!is.null(out_path)) {
      dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
      grDevices::png(filename = out_path, width = width, height = height)
      bfs <- tryCatch(render_snaps(), finally = grDevices::dev.off())
    } else {
      bfs <- render_snaps()
    }
  } else if (mode == "locator") {
    if (!interactive()) {
      stop("mode = 'locator' requires an interactive R session; use mode = 'html' or 'snapshot' non-interactively",
           call. = FALSE)
    }
    bfs <- draw_bfs_panel(source[1L])
    message("Click any vertex to inspect its BFS tree; right-click or press Escape to stop.")
    repeat {
      pt <- graphics::locator(1L)
      if (is.null(pt)) break
      d2 <- (pos[, 1L] - pt$x)^2 + (pos[, 2L] - pt$y)^2
      nearest <- which.min(d2)
      bfs <- draw_bfs_panel(nearest)
    }
  }

  structure(
    list(
      bfs = bfs,
      pos = pos,
      mode = mode,
      output = out_path,
      max_dist = max_dist,
      directed = directed
    ),
    class = "rgraphtool_interactive_bfs"
  )
}

print.rgraphtool_interactive_bfs <- function(x, ...) {
  cat("rgraphtool interactive BFS visualization\n")
  cat(sprintf("  mode: %s%s\n", x$mode,
              if (!is.null(x$output)) paste0(" -> ", x$output) else ""))
  print(x$bfs)
  invisible(x)
}
