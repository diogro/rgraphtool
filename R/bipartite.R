resolve_bipartite_input <- function(graph, type = NULL, dim_labels = NULL,
                                    weight_mode = "auto", weight_name = "weight") {
  if (is.matrix(graph) || inherits(graph, "Matrix") || inherits(graph, "igraph")) {
    graph <- bipartite_graph(graph, type = type, dim_labels = dim_labels,
                             weight_mode = weight_mode, weight_name = weight_name)
    t_int <- validate_bipartite_type(graph, NULL)
  } else {
    t_int <- validate_bipartite_type(graph, type)
    check_bipartite_edges(graph, t_int)
  }
  dl <- get_dim_labels(graph, dim_labels)
  list(graph = graph, type = t_int, dim_labels = dl)
}

resolve_bipartite_weights <- function(graph, weights = NULL, types = NULL,
                                      params = NULL, microcanonical = FALSE) {
  if (is.null(weights)) {
    return(list(
      W = matrix(0, nrow = 0L, ncol = 0L),
      types = character(0),
      params = list(),
      cpp_params = list(),
      weighted = FALSE
    ))
  }
  W <- resolve_weight_matrix(graph, weights)
  p_arg <- if (is.null(params)) {
    if (isTRUE(microcanonical)) "microcanonical" else "bayesian"
  } else {
    params
  }
  spec <- resolve_covariate_specs(W, types, params = p_arg)
  list(
    W = W,
    types = spec$types,
    params = spec$params,
    cpp_params = spec$cpp_params,
    weighted = TRUE
  )
}

canonicalize_bipartite_blocks <- function(blocks, type) {
  n <- length(type)
  if (is.null(blocks)) {
    b0 <- if (any(type == 0L)) 1L else 0L
    out <- integer(n)
    out[type == 0L] <- 1L
    out[type == 1L] <- b0 + 1L
    return(out)
  }
  if (!is.numeric(blocks) || length(blocks) != n || anyNA(blocks) ||
      any(!is.finite(blocks)) || any(blocks < 1) || any(blocks != floor(blocks))) {
    stop("blocks must be a positive integer vector of length graph_num_vertices(graph)",
         call. = FALSE)
  }
  blocks <- as.integer(blocks)
  u0 <- sort(unique(blocks[type == 0L]))
  u1 <- sort(unique(blocks[type == 1L]))
  if (length(intersect(u0, u1)) > 0L) {
    stop("blocks must not assign vertices of different bipartite types to the same block",
         call. = FALSE)
  }
  out <- integer(n)
  if (length(u0)) out[type == 0L] <- match(blocks[type == 0L], u0)
  if (length(u1)) out[type == 1L] <- length(u0) + match(blocks[type == 1L], u1)
  out
}

new_bipartite_state <- function(graph, type, blocks, degree_corrected,
                                degree_prior, edge_prior, wspec, dim_labels,
                                fit = NULL) {
  b0 <- length(unique(blocks[type == 0L]))
  b1 <- length(unique(blocks[type == 1L]))
  B_vec <- stats::setNames(c(b0, b1), dim_labels)
  structure(list(
    blocks = as.integer(blocks),
    type = as.integer(type),
    graph = graph,
    dim_labels = dim_labels,
    degree_corrected = isTRUE(degree_corrected),
    degree_prior = degree_prior,
    edge_prior = edge_prior,
    weights = wspec$W,
    types = wspec$types,
    params = wspec$params,
    cpp_params = wspec$cpp_params,
    weighted = isTRUE(wspec$weighted),
    num_vertices = as.integer(graph_num_vertices(graph)),
    num_edges = as.integer(graph_num_edges(graph)),
    B = B_vec,
    fit = fit
  ), class = "rgraphtool_bipartite_block_state")
}

validate_bipartite_state <- function(state) {
  if (!inherits(state, "rgraphtool_bipartite_block_state")) {
    stop("state must be created by bipartite_block_state() or minimize_bipartite_blockmodel_dl()",
         call. = FALSE)
  }
  if (is.null(state$graph) ||
      graph_num_vertices(state$graph) != state$num_vertices ||
      graph_num_edges(state$graph) != state$num_edges) {
    stop("the graph has changed since the bipartite block state was created", call. = FALSE)
  }
}

#' @noRd
#' @export
bipartite_block_state <- function(graph, type = NULL, blocks = NULL,
                                  degree_corrected = TRUE,
                                  degree_prior = c("distributed", "uniform"),
                                  edge_prior = c("bipartite", "unipartite"),
                                  weights = NULL, types = NULL, params = NULL,
                                  microcanonical = FALSE, dim_labels = NULL) {
  degree_prior <- match.arg(degree_prior)
  edge_prior <- match.arg(edge_prior)
  if (!is.logical(degree_corrected) || length(degree_corrected) != 1L || is.na(degree_corrected)) {
    stop("degree_corrected must be TRUE or FALSE", call. = FALSE)
  }
  inp <- resolve_bipartite_input(graph, type = type, dim_labels = dim_labels)
  b <- canonicalize_bipartite_blocks(blocks, inp$type)
  wspec <- resolve_bipartite_weights(inp$graph, weights = weights, types = types,
                                     params = params, microcanonical = microcanonical)
  new_bipartite_state(inp$graph, inp$type, b, degree_corrected, degree_prior,
                      edge_prior, wspec, inp$dim_labels)
}

#' @noRd
#' @export
bipartite_description_length <- function(state) {
  if (inherits(state, "rgraphtool_nested_bipartite_block_state")) {
    return(nested_bipartite_description_length(state))
  }
  validate_bipartite_state(state)
  cpp_bipartite_entropy(
    graph_pointer(state$graph),
    state$type,
    state$blocks,
    state$degree_corrected,
    identical(state$degree_prior, "distributed"),
    identical(state$edge_prior, "bipartite"),
    state$weights,
    state$types,
    state$cpp_params,
    q_cache_max()
  )
}

#' @noRd
#' @export
minimize_bipartite_blockmodel_dl <- function(graph, type = NULL,
                                             degree_corrected = TRUE,
                                             degree_prior = c("distributed", "uniform"),
                                             edge_prior = c("bipartite", "unipartite"),
                                             weights = NULL, types = NULL, params = NULL,
                                             microcanonical = FALSE,
                                             B_row_range = NULL, B_col_range = NULL,
                                             nruns = 1L, cores = NULL, refine = TRUE,
                                             control = list(),
                                             init_blocks = NULL,
                                             dim_labels = NULL,
                                             weight_mode = c("auto", "binary", "multigraph", "weight")) {
  degree_prior <- match.arg(degree_prior)
  edge_prior <- match.arg(edge_prior)
  weight_mode <- match.arg(weight_mode)
  if (!is.logical(degree_corrected) || length(degree_corrected) != 1L || is.na(degree_corrected)) {
    stop("degree_corrected must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.numeric(nruns) || length(nruns) != 1L || is.na(nruns) || nruns < 1 || nruns != floor(nruns)) {
    stop("nruns must be a positive integer", call. = FALSE)
  }
  nruns <- as.integer(nruns)

  is_mat <- is.matrix(graph) || inherits(graph, "Matrix")
  inp <- resolve_bipartite_input(graph, type = type, dim_labels = dim_labels,
                                 weight_mode = weight_mode)
  if (is_mat && is.null(weights) && "weight" %in% graph_edge_properties(inp$graph)) {
    weights <- "weight"
    if (is.null(types)) types <- "real-normal"
  }
  wspec <- resolve_bipartite_weights(inp$graph, weights = weights, types = types,
                                     params = params, microcanonical = microcanonical)

  n0 <- sum(inp$type == 0L)
  n1 <- sum(inp$type == 1L)
  parse_range <- function(rng, n_side, nm) {
    if (is.null(rng)) return(c(if (n_side > 0L) 1L else 0L, n_side))
    if (!is.numeric(rng) || length(rng) != 2L || anyNA(rng) || rng[1L] < 1 || rng[2L] < rng[1L]) {
      stop(sprintf("%s must be a length-2 integer range c(min, max)", nm), call. = FALSE)
    }
    as.integer(c(max(1L, floor(rng[1L])), min(n_side, floor(rng[2L]))))
  }
  r0 <- parse_range(B_row_range, n0, "B_row_range")
  r1 <- parse_range(B_col_range, n1, "B_col_range")

  r_shrink <- if (!is.null(control$r)) as.numeric(control$r) else 2.0
  n_cand <- if (!is.null(control$n_cand)) as.integer(control$n_cand) else 12L
  init_b <- if (!is.null(init_blocks)) canonicalize_bipartite_blocks(init_blocks, inp$type) else NULL

  best_res <- NULL
  .starts <- parallel_starts(nruns, function(run) {
    res <- cpp_bipartite_minimize(
      graph_pointer(inp$graph),
      inp$type,
      degree_corrected,
      identical(degree_prior, "distributed"),
      identical(edge_prior, "bipartite"),
      wspec$W,
      wspec$types,
      wspec$cpp_params,
      r0[1L], r0[2L],
      r1[1L], r1[2L],
      r_shrink,
      n_cand,
      isTRUE(refine),
      q_cache_max(),
      init_b
    )
    res
  }, cores = cores, graph = graph)
  for (run in seq_len(nruns)) {
    res <- .starts[[run]]
    if (is.null(best_res) || res$description_length < best_res$description_length) {
      best_res <- res
    }
  }

  new_bipartite_state(
    inp$graph, inp$type, best_res$blocks, degree_corrected, degree_prior,
    edge_prior, wspec, inp$dim_labels,
    fit = list(description_length = best_res$description_length, nruns = nruns)
  )
}

# ---------------------------------------------------------------------------
# Nested (hierarchical) bipartite block state

new_nested_bipartite_state <- function(graph, type, levels, level_entropies,
                                       degree_corrected, degree_prior,
                                       edge_prior, wspec, dim_labels) {
  b0 <- if (length(levels)) length(unique(levels[[1L]][type == 0L])) else 0L
  b1 <- if (length(levels)) length(unique(levels[[1L]][type == 1L])) else 0L
  structure(list(
    levels = levels,
    blocks = if (length(levels)) levels[[1L]] else integer(0),
    type = as.integer(type),
    graph = graph,
    dim_labels = dim_labels,
    degree_corrected = isTRUE(degree_corrected),
    degree_prior = degree_prior,
    edge_prior = edge_prior,
    weights = wspec$W,
    types = wspec$types,
    params = wspec$params,
    cpp_params = wspec$cpp_params,
    weighted = isTRUE(wspec$weighted),
    level_entropies = as.numeric(level_entropies),
    description_length = sum(level_entropies),
    num_vertices = as.integer(graph_num_vertices(graph)),
    num_edges = as.integer(graph_num_edges(graph)),
    B = stats::setNames(c(b0, b1), dim_labels)
  ), class = "rgraphtool_nested_bipartite_block_state")
}

validate_nested_bipartite_state <- function(state) {
  if (!inherits(state, "rgraphtool_nested_bipartite_block_state")) {
    stop("state must be created by nested_bipartite_block_state() or minimize_nested_bipartite_blockmodel_dl()",
         call. = FALSE)
  }
  if (is.null(state$graph) ||
      graph_num_vertices(state$graph) != state$num_vertices ||
      graph_num_edges(state$graph) != state$num_edges) {
    stop("the graph has changed since the nested bipartite block state was created", call. = FALSE)
  }
}

#' @noRd
#' @export
nested_bipartite_block_state <- function(graph, type = NULL, levels = NULL,
                                         degree_corrected = TRUE,
                                         degree_prior = c("distributed", "uniform"),
                                         edge_prior = c("bipartite", "unipartite"),
                                         weights = NULL, types = NULL, params = NULL,
                                         microcanonical = FALSE, dim_labels = NULL) {
  degree_prior <- match.arg(degree_prior)
  edge_prior <- match.arg(edge_prior)
  inp <- resolve_bipartite_input(graph, type = type, dim_labels = dim_labels)
  wspec <- resolve_bipartite_weights(inp$graph, weights = weights, types = types,
                                     params = params, microcanonical = microcanonical)
  if (is.null(levels)) {
    levels <- list(canonicalize_bipartite_blocks(NULL, inp$type))
  } else {
    if (!is.list(levels) || !length(levels)) {
      stop("levels must be a non-empty list of integer block vectors", call. = FALSE)
    }
    levels <- lapply(levels, as.integer)
  }
  lev_S <- cpp_nested_bipartite_entropy(
    graph_pointer(inp$graph),
    inp$type,
    levels,
    degree_corrected,
    identical(degree_prior, "distributed"),
    identical(edge_prior, "bipartite"),
    wspec$W,
    wspec$types,
    wspec$cpp_params,
    q_cache_max()
  )
  new_nested_bipartite_state(inp$graph, inp$type, levels, lev_S,
                             degree_corrected, degree_prior, edge_prior,
                             wspec, inp$dim_labels)
}

#' @noRd
#' @export
nested_bipartite_description_length <- function(state) {
  validate_nested_bipartite_state(state)
  sum(cpp_nested_bipartite_entropy(
    graph_pointer(state$graph),
    state$type,
    state$levels,
    state$degree_corrected,
    identical(state$degree_prior, "distributed"),
    identical(state$edge_prior, "bipartite"),
    state$weights,
    state$types,
    state$cpp_params,
    q_cache_max()
  ))
}

#' @noRd
#' @export
minimize_nested_bipartite_blockmodel_dl <- function(graph, type = NULL,
                                                    degree_corrected = TRUE,
                                                    degree_prior = c("distributed", "uniform"),
                                                    edge_prior = c("bipartite", "unipartite"),
                                                    weights = NULL, types = NULL, params = NULL,
                                                    microcanonical = FALSE,
                                                    nruns = 1L, cores = NULL,
                                                    control = list(),
                                                    dim_labels = NULL,
                                                    weight_mode = c("auto", "binary", "multigraph", "weight")) {
  degree_prior <- match.arg(degree_prior)
  edge_prior <- match.arg(edge_prior)
  weight_mode <- match.arg(weight_mode)
  nruns <- as.integer(nruns)
  if (is.na(nruns) || nruns < 1L) stop("nruns must be a positive integer", call. = FALSE)

  is_mat <- is.matrix(graph) || inherits(graph, "Matrix")
  inp <- resolve_bipartite_input(graph, type = type, dim_labels = dim_labels,
                                 weight_mode = weight_mode)
  if (is_mat && is.null(weights) && "weight" %in% graph_edge_properties(inp$graph)) {
    weights <- "weight"
    if (is.null(types)) types <- "real-normal"
  }
  wspec <- resolve_bipartite_weights(inp$graph, weights = weights, types = types,
                                     params = params, microcanonical = microcanonical)

  r_shrink <- if (!is.null(control$r)) as.numeric(control$r) else 2.0
  n_cand <- if (!is.null(control$n_cand)) as.integer(control$n_cand) else 12L

  best_res <- NULL
  .starts <- parallel_starts(nruns, function(run) {
    res <- cpp_nested_bipartite_minimize(
      graph_pointer(inp$graph),
      inp$type,
      degree_corrected,
      identical(degree_prior, "distributed"),
      identical(edge_prior, "bipartite"),
      wspec$W,
      wspec$types,
      wspec$cpp_params,
      r_shrink,
      n_cand,
      q_cache_max()
    )
    res
  }, cores = cores, graph = graph)
  for (run in seq_len(nruns)) {
    res <- .starts[[run]]
    if (is.null(best_res) || res$description_length < best_res$description_length) {
      best_res <- res
    }
  }

  new_nested_bipartite_state(
    inp$graph, inp$type, best_res$levels, best_res$level_entropies,
    degree_corrected, degree_prior, edge_prior, wspec, inp$dim_labels
  )
}

# ---------------------------------------------------------------------------
# Inspection helpers: memberships, rectangular block edge counts & affinity

#' @noRd
#' @export
bipartite_membership <- function(state, level = 1L) {
  if (inherits(state, "rgraphtool_nested_bipartite_block_state")) {
    validate_nested_bipartite_state(state)
    b <- state$levels[[1L]]
    if (level > 1L) {
      if (level > length(state$levels)) {
        stop("level exceeds hierarchy depth", call. = FALSE)
      }
      for (l in 2L:level) b <- state$levels[[l]][b]
    }
  } else {
    validate_bipartite_state(state)
    b <- state$blocks
  }

  rows <- which(state$type == 0L)
  cols <- which(state$type == 1L)
  b_row <- match(b[rows], sort(unique(b[rows])))
  b_col <- match(b[cols], sort(unique(b[cols])))

  if ("name" %in% graph_vertex_properties(state$graph)) {
    nms <- graph_vertex_property(state$graph, "name")
    if (is.character(nms) && length(nms) == length(state$type)) {
      names(b_row) <- nms[rows]
      names(b_col) <- nms[cols]
    }
  }
  stats::setNames(list(b_row, b_col), state$dim_labels)
}

#' @noRd
#' @export
bipartite_edge_counts <- function(state, level = 1L) {
  mem <- bipartite_membership(state, level = level)
  b_row <- mem[[1L]]
  b_col <- mem[[2L]]
  B0 <- if (length(b_row)) max(b_row) else 0L
  B1 <- if (length(b_col)) max(b_col) else 0L
  mat <- matrix(0, nrow = B0, ncol = B1)

  rows <- which(state$type == 0L)
  cols <- which(state$type == 1L)
  row_pos <- integer(length(state$type))
  col_pos <- integer(length(state$type))
  row_pos[rows] <- seq_along(rows)
  col_pos[cols] <- seq_along(cols)

  el <- graph_edges(state$graph)
  if (nrow(el) > 0L) {
    for (i in seq_len(nrow(el))) {
      u <- el[i, 1L]
      v <- el[i, 2L]
      if (state$type[u] == 0L) {
        r <- b_row[row_pos[u]]
        s <- b_col[col_pos[v]]
      } else {
        r <- b_row[row_pos[v]]
        s <- b_col[col_pos[u]]
      }
      mat[r, s] <- mat[r, s] + 1
    }
  }
  dl <- state$dim_labels
  dimnames(mat) <- stats::setNames(
    list(paste0(dl[1L], "_", seq_len(B0)), paste0(dl[2L], "_", seq_len(B1))),
    dl
  )
  mat
}

#' @noRd
#' @export
bipartite_block_matrix <- function(state, level = 1L,
                                   normalize = c("density", "counts", "mean_weight")) {
  normalize <- match.arg(normalize)
  counts <- bipartite_edge_counts(state, level = level)
  if (normalize == "counts") return(counts)

  mem <- bipartite_membership(state, level = level)
  n_row <- as.numeric(tabulate(mem[[1L]], nbins = nrow(counts)))
  n_col <- as.numeric(tabulate(mem[[2L]], nbins = ncol(counts)))

  if (normalize == "density") {
    denom <- outer(n_row, n_col)
    out <- ifelse(denom > 0, counts / denom, 0)
    dimnames(out) <- dimnames(counts)
    return(out)
  }

  # mean_weight (first covariate column)
  if (!isTRUE(state$weighted) || ncol(state$weights) < 1L) {
    stop("normalize = 'mean_weight' requires a weighted bipartite block state",
         call. = FALSE)
  }
  rows <- which(state$type == 0L)
  cols <- which(state$type == 1L)
  row_pos <- integer(length(state$type))
  col_pos <- integer(length(state$type))
  row_pos[rows] <- seq_along(rows)
  col_pos[cols] <- seq_along(cols)

  sum_w <- matrix(0, nrow = nrow(counts), ncol = ncol(counts))
  el <- graph_edges(state$graph)
  w <- state$weights[, 1L]
  if (nrow(el) > 0L) {
    for (i in seq_len(nrow(el))) {
      u <- el[i, 1L]
      v <- el[i, 2L]
      if (state$type[u] == 0L) {
        r <- mem[[1L]][row_pos[u]]
        s <- mem[[2L]][col_pos[v]]
      } else {
        r <- mem[[1L]][row_pos[v]]
        s <- mem[[2L]][col_pos[u]]
      }
      sum_w[r, s] <- sum_w[r, s] + w[i]
    }
  }
  out <- ifelse(counts > 0, sum_w / counts, NA_real_)
  dimnames(out) <- dimnames(counts)
  out
}

#' @noRd
#' @export
compare_bipartite_models <- function(graph, type = NULL,
                                     weights = NULL, types = NULL, params = NULL,
                                     microcanonical = FALSE,
                                     edge_prior = c("bipartite", "unipartite"),
                                     nested = TRUE, nruns = 3L,
                                     dim_labels = NULL) {
  edge_prior <- match.arg(edge_prior)
  inp <- resolve_bipartite_input(graph, type = type, dim_labels = dim_labels)

  fits <- list()
  rows <- list()
  add_model <- function(label, st) {
    dl <- bipartite_description_length(st)
    fits[[label]] <<- st
    rows[[label]] <<- data.frame(
      model = label,
      nested = inherits(st, "rgraphtool_nested_bipartite_block_state"),
      degree_corrected = st$degree_corrected,
      B_row = unname(st$B[1L]),
      B_col = unname(st$B[2L]),
      description_length = dl,
      stringsAsFactors = FALSE
    )
  }

  add_model("flat_dc", minimize_bipartite_blockmodel_dl(
    inp$graph, type = inp$type, degree_corrected = TRUE, edge_prior = edge_prior,
    weights = weights, types = types, params = params, microcanonical = microcanonical,
    nruns = nruns, dim_labels = inp$dim_labels
  ))
  add_model("flat_ndc", minimize_bipartite_blockmodel_dl(
    inp$graph, type = inp$type, degree_corrected = FALSE, edge_prior = edge_prior,
    weights = weights, types = types, params = params, microcanonical = microcanonical,
    nruns = nruns, dim_labels = inp$dim_labels
  ))
  if (isTRUE(nested)) {
    add_model("nested_dc", minimize_nested_bipartite_blockmodel_dl(
      inp$graph, type = inp$type, degree_corrected = TRUE, edge_prior = edge_prior,
      weights = weights, types = types, params = params, microcanonical = microcanonical,
      nruns = nruns, dim_labels = inp$dim_labels
    ))
    add_model("nested_ndc", minimize_nested_bipartite_blockmodel_dl(
      inp$graph, type = inp$type, degree_corrected = FALSE, edge_prior = edge_prior,
      weights = weights, types = types, params = params, microcanonical = microcanonical,
      nruns = nruns, dim_labels = inp$dim_labels
    ))
  }

  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL
  tab$delta_dl <- tab$description_length - min(tab$description_length)
  tab <- tab[order(tab$description_length), , drop = FALSE]
  rownames(tab) <- NULL
  structure(list(
    comparison = tab,
    best_model = tab$model[1L],
    fits = fits
  ), class = "rgraphtool_bipartite_model_comparison")
}

# ---------------------------------------------------------------------------
# Print methods

#' @noRd
#' @export
print.rgraphtool_bipartite_block_state <- function(x, ...) {
  dl <- if (!is.null(x$fit$description_length)) {
    x$fit$description_length
  } else {
    bipartite_description_length(x)
  }
  n0 <- sum(x$type == 0L)
  n1 <- sum(x$type == 1L)
  cat(sprintf("<bipartite block state (%s, %s, edge_prior = %s)>\n",
              if (x$degree_corrected) "degree-corrected" else "non-degree-corrected",
              if (x$weighted) paste("weighted:", paste(x$types, collapse = ", ")) else "unweighted",
              x$edge_prior))
  cat(sprintf("  %s: %d vertices -> %d blocks\n", x$dim_labels[1L], n0, unname(x$B[1L])))
  cat(sprintf("  %s: %d vertices -> %d blocks\n", x$dim_labels[2L], n1, unname(x$B[2L])))
  cat(sprintf("  edges: %d | description length: %.4f nats\n", x$num_edges, dl))
  invisible(x)
}

#' @noRd
#' @export
print.rgraphtool_nested_bipartite_block_state <- function(x, ...) {
  n0 <- sum(x$type == 0L)
  n1 <- sum(x$type == 1L)
  cat(sprintf("<nested bipartite block state (%s, %d levels, edge_prior = %s)>\n",
              if (x$degree_corrected) "degree-corrected" else "non-degree-corrected",
              length(x$levels), x$edge_prior))
  cat(sprintf("  %s: %d vertices -> %d base blocks\n", x$dim_labels[1L], n0, unname(x$B[1L])))
  cat(sprintf("  %s: %d vertices -> %d base blocks\n", x$dim_labels[2L], n1, unname(x$B[2L])))
  cat(sprintf("  edges: %d | description length: %.4f nats\n",
              x$num_edges, x$description_length))
  invisible(x)
}

#' @noRd
#' @export
print.rgraphtool_bipartite_model_comparison <- function(x, ...) {
  cat("<bipartite block model comparison>\n")
  print(x$comparison, row.names = FALSE)
  cat(sprintf("Best model: %s\n", x$best_model))
  invisible(x)
}
