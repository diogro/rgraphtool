resolve_ranked_u <- function(u, b, unique, distance, M) {
  B <- max(b, 0L)
  if (is.null(u)) {
    u_vec <- seq_len(B) - 1
  } else {
    if (!is.numeric(u) || anyNA(u) || any(!is.finite(u)) || any(u != floor(u))) {
      stop("u must be a finite integer vector of group ranks", call. = FALSE)
    }
    if (length(u) < B) {
      stop(sprintf("u must have length at least B = %d", B), call. = FALSE)
    }
    u_vec <- as.numeric(u[seq_len(B)])
  }
  if (isTRUE(unique) && length(unique(u_vec)) < length(u_vec)) {
    stop("u must contain distinct ranks when unique = TRUE", call. = FALSE)
  }
  if (isTRUE(distance) && isTRUE(unique)) {
    u_vec <- seq_len(B) - 1
  }
  u_vec
}

new_ranked_state <- function(graph, blocks, u, distance, unique, M,
                             degree_corrected, degree_prior, fit = NULL) {
  structure(list(
    blocks = as.integer(blocks),
    u = as.numeric(u),
    graph = graph,
    distance = isTRUE(distance),
    unique = isTRUE(unique),
    M = as.integer(M),
    degree_corrected = isTRUE(degree_corrected),
    degree_prior = degree_prior,
    num_vertices = graph_num_vertices(graph),
    num_edges = graph_num_edges(graph),
    B = length(unique(blocks)),
    fit = fit
  ), class = "graphtool_ranked_block_state")
}

validate_ranked_state <- function(state) {
  if (!inherits(state, "graphtool_ranked_block_state")) {
    stop("state must be created by ranked_block_state() or minimize_ranked_blockmodel_dl()",
         call. = FALSE)
  }
  if (graph_num_edges(state$graph) != state$num_edges ||
      graph_num_vertices(state$graph) != state$num_vertices) {
    stop("the graph has changed since the ranked block state was created",
         call. = FALSE)
  }
}

ranked_block_state <- function(graph, blocks = NULL, u = NULL,
                               distance = FALSE, unique = FALSE, M = NULL,
                               degree_corrected = TRUE,
                               degree_prior = c("distributed", "uniform")) {
  graph_pointer(graph)
  if (!isTRUE(graph$directed)) {
    stop("ranked_block_state() requires a directed graph", call. = FALSE)
  }
  degree_prior <- match.arg(degree_prior)
  if (!is.logical(distance) || length(distance) != 1L || is.na(distance)) {
    stop("distance must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(unique) || length(unique) != 1L || is.na(unique)) {
    stop("unique must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(degree_corrected) || length(degree_corrected) != 1L ||
      is.na(degree_corrected)) {
    stop("degree_corrected must be TRUE or FALSE", call. = FALSE)
  }
  n <- graph_num_vertices(graph)
  M_val <- if (is.null(M)) max(1L, as.integer(n)) else as.integer(M)
  if (is.na(M_val) || M_val < 1L) stop("M must be a positive integer", call. = FALSE)

  if (is.null(blocks)) {
    blocks <- rep.int(1L, n)
  } else {
    blocks <- validate_label_vector(blocks, "blocks")
    if (length(blocks) != n || anyNA(blocks)) {
      stop("blocks must be an integer vector of length graph_num_vertices(graph)",
           call. = FALSE)
    }
  }
  u_vec <- resolve_ranked_u(u, blocks, unique, distance, M_val)
  new_ranked_state(graph, blocks, u_vec, distance, unique, M_val,
                   degree_corrected, degree_prior)
}

ranked_description_length <- function(state, label_dl = TRUE, degree_prior = NULL) {
  validate_ranked_state(state)
  if (is.null(degree_prior)) degree_prior <- state$degree_prior
  degree_prior <- match.arg(degree_prior, c("distributed", "uniform"))
  cpp_ranked_entropy(graph_pointer(state$graph), state$blocks, state$u,
                     state$distance, state$unique, as.numeric(state$M),
                     state$degree_corrected, degree_prior == "distributed",
                     isTRUE(label_dl), q_cache_max())
}

minimize_ranked_blockmodel_dl <- function(graph, distance = FALSE, unique = FALSE,
                                          M = NULL, degree_corrected = TRUE,
                                          degree_prior = c("distributed", "uniform"),
                                          B_min = 1L, B_max = Inf,
                                          refine = TRUE, epsilon = 1e-3) {
  pointer <- graph_pointer(graph)
  if (!isTRUE(graph$directed)) {
    stop("minimize_ranked_blockmodel_dl() requires a directed graph", call. = FALSE)
  }
  degree_prior <- match.arg(degree_prior)
  if (!is.logical(distance) || length(distance) != 1L || is.na(distance)) {
    stop("distance must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(unique) || length(unique) != 1L || is.na(unique)) {
    stop("unique must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(degree_corrected) || length(degree_corrected) != 1L ||
      is.na(degree_corrected)) {
    stop("degree_corrected must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.numeric(B_min) || length(B_min) != 1L || is.na(B_min) || B_min < 1 ||
      !is.numeric(B_max) || length(B_max) != 1L || is.na(B_max) || B_max < B_min) {
    stop("B_min and B_max must satisfy 1 <= B_min <= B_max", call. = FALSE)
  }
  n <- graph_num_vertices(graph)
  M_val <- if (is.null(M)) max(1L, as.integer(n)) else as.integer(M)
  if (is.na(M_val) || M_val < 1L) stop("M must be a positive integer", call. = FALSE)

  res <- cpp_ranked_minimize(pointer, distance, unique, as.numeric(M_val),
                             degree_corrected, degree_prior == "distributed",
                             as.integer(floor(B_min)),
                             if (is.finite(B_max)) floor(B_max) else Inf,
                             isTRUE(refine), as.numeric(epsilon), q_cache_max())
  new_ranked_state(graph, res$blocks, res$u, distance, unique, M_val,
                   degree_corrected, degree_prior,
                   fit = list(description_length = res$entropy))
}

print.graphtool_ranked_block_state <- function(x, ...) {
  dl <- ranked_description_length(x)
  ec <- ranked_edge_counts(x)
  cat(sprintf("<ranked block state: %s, B = %d, %d upstream / %d lateral / %d downstream, DL = %.4f nats>\n",
              if (x$degree_corrected) "degree-corrected" else "non-degree-corrected",
              x$B, ec$upstream, ec$lateral, ec$downstream, dl))
  invisible(x)
}

# ---------------------------------------------------------------------------
# Nested Ranked SBM (NestedBlockState(base_state=RankedBlockState))

new_nested_ranked_state <- function(graph, levels, u, distance, unique, M,
                                    degree_corrected, degree_prior, fit = NULL) {
  structure(list(
    levels = levels,
    u = as.numeric(u),
    graph = graph,
    distance = isTRUE(distance),
    unique = isTRUE(unique),
    M = as.integer(M),
    degree_corrected = isTRUE(degree_corrected),
    degree_prior = degree_prior,
    num_vertices = graph_num_vertices(graph),
    num_edges = graph_num_edges(graph),
    B = vapply(levels, function(b) length(unique(b)), integer(1)),
    fit = fit
  ), class = "graphtool_nested_ranked_block_state")
}

validate_nested_ranked_state <- function(state) {
  if (!inherits(state, "graphtool_nested_ranked_block_state")) {
    stop("state must be created by nested_ranked_block_state() or minimize_nested_ranked_blockmodel_dl()",
         call. = FALSE)
  }
  if (graph_num_edges(state$graph) != state$num_edges ||
      graph_num_vertices(state$graph) != state$num_vertices) {
    stop("the graph has changed since the nested ranked block state was created",
         call. = FALSE)
  }
}

nested_ranked_block_state <- function(graph, levels = NULL, u = NULL,
                                      distance = FALSE, unique = FALSE, M = NULL,
                                      degree_corrected = TRUE,
                                      degree_prior = c("distributed", "uniform")) {
  graph_pointer(graph)
  if (!isTRUE(graph$directed)) {
    stop("nested_ranked_block_state() requires a directed graph", call. = FALSE)
  }
  degree_prior <- match.arg(degree_prior)
  if (!is.logical(distance) || length(distance) != 1L || is.na(distance)) {
    stop("distance must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(unique) || length(unique) != 1L || is.na(unique)) {
    stop("unique must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(degree_corrected) || length(degree_corrected) != 1L ||
      is.na(degree_corrected)) {
    stop("degree_corrected must be TRUE or FALSE", call. = FALSE)
  }
  n <- graph_num_vertices(graph)
  if (n < 1L) stop("nested ranked states need at least one vertex", call. = FALSE)
  M_val <- if (is.null(M)) max(1L, as.integer(n)) else as.integer(M)
  if (is.na(M_val) || M_val < 1L) stop("M must be a positive integer", call. = FALSE)

  if (is.null(levels)) {
    levels <- list(rep.int(1L, n), 1L)
  } else {
    if (!is.list(levels) || length(levels) < 2L) {
      stop("levels must be a list of at least 2 positive integer label vectors",
           call. = FALSE)
    }
    cur_n <- n
    levels <- lapply(seq_along(levels), function(l) {
      v <- validate_label_vector(levels[[l]], sprintf("levels[[%d]]", l))
      if (length(v) != cur_n || anyNA(v)) {
        stop(sprintf("levels[[%d]] must have %d non-missing positive integer labels",
                     l, cur_n), call. = FALSE)
      }
      cur_n <<- length(unique(v))
      v
    })
  }
  u_vec <- resolve_ranked_u(u, levels[[1L]], unique, distance, M_val)
  new_nested_ranked_state(graph, levels, u_vec, distance, unique, M_val,
                          degree_corrected, degree_prior)
}

nested_ranked_level_entropy <- function(state, degree_prior = NULL) {
  validate_nested_ranked_state(state)
  if (is.null(degree_prior)) degree_prior <- state$degree_prior
  degree_prior <- match.arg(degree_prior, c("distributed", "uniform"))
  cpp_nested_ranked_entropy(graph_pointer(state$graph), state$levels, state$u,
                            state$distance, state$unique, as.numeric(state$M),
                            state$degree_corrected, degree_prior == "distributed",
                            q_cache_max())
}

nested_ranked_description_length <- function(state, degree_prior = NULL) {
  sum(nested_ranked_level_entropy(state, degree_prior = degree_prior))
}

minimize_nested_ranked_blockmodel_dl <- function(graph, distance = FALSE,
                                                 unique = FALSE, M = NULL,
                                                 degree_corrected = TRUE,
                                                 degree_prior = c("distributed", "uniform"),
                                                 B_min = 1L, B_max = Inf,
                                                 refine = TRUE, epsilon = 1e-3) {
  pointer <- graph_pointer(graph)
  if (!isTRUE(graph$directed)) {
    stop("minimize_nested_ranked_blockmodel_dl() requires a directed graph",
         call. = FALSE)
  }
  degree_prior <- match.arg(degree_prior)
  if (graph_num_vertices(graph) < 1L) {
    stop("nested ranked states need at least one vertex", call. = FALSE)
  }
  if (!is.logical(distance) || length(distance) != 1L || is.na(distance)) {
    stop("distance must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(unique) || length(unique) != 1L || is.na(unique)) {
    stop("unique must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(degree_corrected) || length(degree_corrected) != 1L ||
      is.na(degree_corrected)) {
    stop("degree_corrected must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.numeric(B_min) || length(B_min) != 1L || is.na(B_min) || B_min < 1 ||
      !is.numeric(B_max) || length(B_max) != 1L || is.na(B_max) || B_max < B_min) {
    stop("B_min and B_max must satisfy 1 <= B_min <= B_max", call. = FALSE)
  }
  n <- graph_num_vertices(graph)
  M_val <- if (is.null(M)) max(1L, as.integer(n)) else as.integer(M)
  if (is.na(M_val) || M_val < 1L) stop("M must be a positive integer", call. = FALSE)

  res <- cpp_nested_ranked_minimize(pointer, distance, unique, as.numeric(M_val),
                                    degree_corrected, degree_prior == "distributed",
                                    as.integer(floor(B_min)),
                                    if (is.finite(B_max)) floor(B_max) else Inf,
                                    isTRUE(refine), as.numeric(epsilon),
                                    q_cache_max())
  new_nested_ranked_state(graph, res$bs, res$u, distance, unique, M_val,
                          degree_corrected, degree_prior,
                          fit = list(description_length = res$entropy,
                                     level_entropy = res$level_entropy))
}

print.graphtool_nested_ranked_block_state <- function(x, ...) {
  dl <- nested_ranked_description_length(x)
  ec <- ranked_edge_counts(x)
  cat(sprintf("<nested ranked block state: %s, %d levels, B = %s, %d up / %d lat / %d down, DL = %.4f nats>\n",
              if (x$degree_corrected) "degree-corrected" else "non-degree-corrected",
              length(x$levels), paste(x$B, collapse = "-"),
              ec$upstream, ec$lateral, ec$downstream, dl))
  invisible(x)
}

# ---------------------------------------------------------------------------
# Inspection helpers for flat and nested ranked states

ranked_base_blocks <- function(state) {
  if (inherits(state, "graphtool_nested_ranked_block_state")) {
    validate_nested_ranked_state(state)
    return(state$levels[[1L]])
  }
  validate_ranked_state(state)
  state$blocks
}

ranked_vertex_order <- function(state, normalize = FALSE) {
  b <- ranked_base_blocks(state)
  u_v <- state$u[b]
  if (isTRUE(normalize)) {
    rng <- range(u_v)
    if (rng[2L] > rng[1L]) {
      return((u_v - rng[1L]) / (rng[2L] - rng[1L]))
    }
    return(rep(0, length(u_v)))
  }
  u_v
}

ranked_edge_directions <- function(state) {
  b <- ranked_base_blocks(state)
  ed <- graph_edges(state$graph)
  if (!nrow(ed)) return(integer(0))
  us <- state$u[b[ed[, 1L]]]
  ut <- state$u[b[ed[, 2L]]]
  ifelse(us < ut, 1L, ifelse(us > ut, -1L, 0L))
}

ranked_edge_counts <- function(state) {
  b <- ranked_base_blocks(state)
  ed <- graph_edges(state$graph)
  if (!nrow(ed)) {
    return(list(upstream = 0L, lateral = 0L, downstream = 0L, internal = 0L))
  }
  bs <- b[ed[, 1L]]
  bt <- b[ed[, 2L]]
  us <- state$u[bs]
  ut <- state$u[bt]
  list(
    upstream = sum(us < ut),
    lateral = sum(us == ut),
    downstream = sum(us > ut),
    internal = sum(bs == bt)
  )
}

# Declarative comparison of ordered (ranked) vs unordered SBMs
compare_ordered_models <- function(graph,
                                   models = c("unordered_flat", "unordered_nested",
                                              "ranked_flat", "ranked_nested"),
                                   degree_corrected = TRUE,
                                   degree_prior = c("distributed", "uniform"),
                                   method = c("map", "posterior"),
                                   evidence_method = c("variational", "importance"),
                                   nruns = 1L, nsamples = 200L, thin = 5L,
                                   burnin = NULL, mle = TRUE,
                                   B_min = 1L, B_max = Inf, refine = TRUE) {
  graph_pointer(graph)
  if (!isTRUE(graph$directed)) {
    stop("compare_ordered_models() requires a directed graph", call. = FALSE)
  }
  degree_prior <- match.arg(degree_prior)
  method <- match.arg(method)
  evidence_method <- match.arg(evidence_method)
  if (!is.logical(degree_corrected) || length(degree_corrected) != 1L ||
      is.na(degree_corrected)) {
    stop("degree_corrected must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.numeric(nruns) || length(nruns) != 1L || is.na(nruns) ||
      nruns < 1 || nruns != floor(nruns)) {
    stop("nruns must be a positive integer", call. = FALSE)
  }
  nruns <- as.integer(nruns)
  valid <- c("unordered_flat", "unordered_nested", "ranked_flat",
             "ranked_nested", "ranked_distance_flat", "ranked_distance_nested")
  if (!is.character(models) || !length(models) || anyNA(models) ||
      any(!models %in% valid)) {
    stop("models must be a non-empty subset of: ", paste(valid, collapse = ", "),
         call. = FALSE)
  }
  models <- unique(models)

  states <- stats::setNames(vector("list", length(models)), models)
  evidences <- if (method == "posterior") {
    stats::setNames(vector("list", length(models)), models)
  } else {
    NULL
  }
  rows <- vector("list", length(models))

  for (i in seq_along(models)) {
    nm <- models[i]
    best_st <- NULL
    best_dl <- Inf
    for (run in seq_len(nruns)) {
      st <- switch(
        nm,
        unordered_flat = minimize_blockmodel_dl(
          graph, degree_corrected = degree_corrected,
          degree_prior = degree_prior, B_min = B_min, B_max = B_max, refine = refine
        ),
        unordered_nested = minimize_nested_blockmodel_dl(
          graph, degree_corrected = degree_corrected,
          degree_prior = degree_prior, B_min = B_min, B_max = B_max, refine = refine
        ),
        ranked_flat = minimize_ranked_blockmodel_dl(
          graph, distance = FALSE, unique = FALSE,
          degree_corrected = degree_corrected, degree_prior = degree_prior,
          B_min = B_min, B_max = B_max, refine = refine
        ),
        ranked_nested = minimize_nested_ranked_blockmodel_dl(
          graph, distance = FALSE, unique = FALSE,
          degree_corrected = degree_corrected, degree_prior = degree_prior,
          B_min = B_min, B_max = B_max, refine = refine
        ),
        ranked_distance_flat = minimize_ranked_blockmodel_dl(
          graph, distance = TRUE, unique = FALSE,
          degree_corrected = degree_corrected, degree_prior = degree_prior,
          B_min = B_min, B_max = B_max, refine = refine
        ),
        ranked_distance_nested = minimize_nested_ranked_blockmodel_dl(
          graph, distance = TRUE, unique = FALSE,
          degree_corrected = degree_corrected, degree_prior = degree_prior,
          B_min = B_min, B_max = B_max, refine = refine
        )
      )
      dl <- state_description_length(st)
      if (dl < best_dl) {
        best_dl <- dl
        best_st <- st
      }
    }
    states[[nm]] <- best_st

    if (method == "posterior") {
      ev <- model_evidence(best_st, nsamples = nsamples, thin = thin,
                           burnin = burnin, mle = mle, method = evidence_method,
                           B_min = B_min, B_max = B_max)
      evidences[[nm]] <- ev
      neg_log_ev <- ev$neg_log_evidence
      post_h <- ev$posterior_entropy
    } else {
      neg_log_ev <- best_dl
      post_h <- NA_real_
    }

    rows[[i]] <- data.frame(
      model = nm,
      B = state_num_groups(best_st),
      description_length = best_dl,
      posterior_entropy = post_h,
      neg_log_evidence = neg_log_ev,
      stringsAsFactors = FALSE
    )
  }

  tab <- do.call(rbind, rows)
  ord <- order(tab$neg_log_evidence, tab$description_length)
  tab <- tab[ord, , drop = FALSE]
  rownames(tab) <- NULL

  scores <- tab$neg_log_evidence
  tab$delta <- scores - min(scores)
  lw <- -tab$delta
  tab$posterior_weight <- exp(lw - logsumexp_vec(lw))

  k <- nrow(tab)
  lbf <- matrix(0, nrow = k, ncol = k, dimnames = list(tab$model, tab$model))
  for (i in seq_len(k)) {
    for (j in seq_len(k)) lbf[i, j] <- scores[j] - scores[i]
  }

  best_model <- tab$model[1L]
  structure(list(
    summary = tab,
    best_model = best_model,
    best_state = states[[best_model]],
    states = states[tab$model],
    evidences = if (is.null(evidences)) NULL else evidences[tab$model],
    log_bayes_factors = lbf,
    degree_corrected = degree_corrected,
    degree_prior = degree_prior,
    method = method,
    evidence_method = if (method == "posterior") evidence_method else NULL
  ), class = "graphtool_ordered_model_comparison")
}

print.graphtool_ordered_model_comparison <- function(x, ...) {
  cat(sprintf("<ordered (ranked) vs. unordered SBM comparison (%s, %s)>\n",
              if (x$degree_corrected) "DC" else "NDC",
              if (x$method == "posterior") paste0("posterior/", x$evidence_method) else "MAP"))
  df <- x$summary
  cols <- if (x$method == "posterior") {
    c("model", "B", "description_length", "posterior_entropy",
      "neg_log_evidence", "delta", "posterior_weight")
  } else {
    c("model", "B", "description_length", "delta", "posterior_weight")
  }
  out <- df[, cols, drop = FALSE]
  for (col in intersect(c("description_length", "posterior_entropy",
                          "neg_log_evidence", "delta"), names(out))) {
    out[[col]] <- round(out[[col]], 3)
  }
  out$posterior_weight <- signif(out$posterior_weight, 3)
  print(out, row.names = FALSE)
  cat(sprintf("Best model: %s\n", x$best_model))
  invisible(x)
}
