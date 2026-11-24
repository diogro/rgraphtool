valid_covariate_types <- function() {
  c("real-exponential", "real-normal", "discrete-geometric",
    "discrete-poisson", "discrete-binomial")
}

resolve_weight_matrix <- function(graph, weights) {
  e <- graph_num_edges(graph)
  if (is.character(weights)) {
    if (!length(weights) || anyNA(weights) || any(!nzchar(weights))) {
      stop("weights character vector must contain non-empty edge property names",
           call. = FALSE)
    }
    cols <- lapply(weights, function(nm) graph_edge_property(graph, nm))
    mat <- do.call(cbind, cols)
    colnames(mat) <- weights
  } else if (is.matrix(weights)) {
    if (!is.numeric(weights) || nrow(weights) != e || ncol(weights) < 1L) {
      stop("weights matrix must be numeric with nrow(weights) == graph_num_edges(graph)",
           call. = FALSE)
    }
    mat <- matrix(as.numeric(weights), nrow = e, ncol = ncol(weights))
    colnames(mat) <- colnames(weights)
  } else if (is.numeric(weights)) {
    if (length(weights) != e) {
      stop("weights vector must have length graph_num_edges(graph)", call. = FALSE)
    }
    mat <- matrix(as.numeric(weights), nrow = e, ncol = 1L)
  } else if (is.list(weights) && length(weights) >= 1L) {
    cols <- lapply(weights, function(item) {
      if (is.character(item) && length(item) == 1L) {
        graph_edge_property(graph, item)
      } else if (is.numeric(item) && length(item) == e) {
        as.numeric(item)
      } else {
        stop("each element of weights list must be an edge property name or numeric vector of length E",
             call. = FALSE)
      }
    })
    mat <- do.call(cbind, cols)
    colnames(mat) <- names(weights)
  } else {
    stop("weights must be an edge property name, numeric vector, E x K matrix, or list",
         call. = FALSE)
  }
  if (anyNA(mat) || any(!is.finite(mat))) {
    stop("weights must contain finite non-missing numeric values", call. = FALSE)
  }
  mat
}

resolve_covariate_specs <- function(weight_mat, types, params = "microcanonical") {
  k <- ncol(weight_mat)
  if (!is.character(types) || length(types) != k || anyNA(types)) {
    stop("types must be a character vector with one entry per covariate", call. = FALSE)
  }
  valid <- valid_covariate_types()
  bad <- setdiff(types, valid)
  if (length(bad)) {
    stop("unknown covariate type(s): ", paste(bad, collapse = ", "),
         ". Valid choices are: ", paste(valid, collapse = ", "), call. = FALSE)
  }

  # Normalize params into a list of length k
  if (is.character(params) && length(params) == 1L) {
    plist <- rep(list(params), k)
  } else if (is.character(params) && length(params) == k) {
    plist <- as.list(params)
  } else if (is.list(params) && k == 1L && !is.null(names(params)) &&
             all(nzchar(names(params)))) {
    plist <- list(params)
  } else if (is.list(params) && length(params) == k) {
    plist <- params
  } else {
    stop("params must be 'microcanonical', 'bayesian', a named parameter list, or a list of length K",
         call. = FALSE)
  }

  cpp_params <- vector("list", k)
  named_params <- vector("list", k)
  for (i in seq_len(k)) {
    w <- weight_mat[, i]
    tp <- types[i]
    pi <- plist[[i]]

    # Validate weight domain
    if (tp == "real-exponential" && any(w < 0)) {
      stop("real-exponential weights must be non-negative", call. = FALSE)
    }
    if (tp %in% c("discrete-geometric", "discrete-poisson", "discrete-binomial")) {
      if (any(w < 0) || any(w != floor(w))) {
        stop(sprintf("%s weights must be non-negative integers", tp), call. = FALSE)
      }
    }

    is_micro <- identical(pi, "microcanonical")
    is_bayes_default <- identical(pi, "bayesian")
    if (!is_micro && !is_bayes_default && !is.list(pi)) {
      stop("each params entry must be 'microcanonical', 'bayesian', or a named list",
           call. = FALSE)
    }

    if (tp %in% c("real-exponential", "discrete-poisson")) {
      if (is_micro) {
        d <- list(alpha = NaN, beta = NaN)
      } else {
        mu <- if (length(w)) mean(w) else 1
        d <- list(alpha = 1, beta = if (mu > 0) mu else 1)
        if (is.list(pi)) {
          unknown <- setdiff(names(pi), names(d))
          if (length(unknown)) stop("unknown parameter(s): ", paste(unknown, collapse = ", "), call. = FALSE)
          for (nm in names(pi)) d[[nm]] <- as.numeric(pi[[nm]])
        }
      }
      cpp_params[[i]] <- c(d$alpha, d$beta)
      named_params[[i]] <- d
    } else if (tp == "discrete-geometric") {
      if (is_micro) {
        d <- list(alpha = NaN, beta = NaN)
      } else {
        d <- list(alpha = 1, beta = 1)
        if (is.list(pi)) {
          unknown <- setdiff(names(pi), names(d))
          if (length(unknown)) stop("unknown parameter(s): ", paste(unknown, collapse = ", "), call. = FALSE)
          for (nm in names(pi)) d[[nm]] <- as.numeric(pi[[nm]])
        }
      }
      cpp_params[[i]] <- c(d$alpha, d$beta)
      named_params[[i]] <- d
    } else if (tp == "discrete-binomial") {
      n_def <- if (length(w)) max(w) else 1
      if (is_micro) {
        d <- list(N = n_def, alpha = NaN, beta = NaN)
      } else if (is.list(pi) && identical(names(pi), "N")) {
        d <- list(N = as.numeric(pi$N), alpha = NaN, beta = NaN)
      } else {
        d <- list(N = n_def, alpha = 1, beta = 1)
        if (is.list(pi)) {
          unknown <- setdiff(names(pi), names(d))
          if (length(unknown)) stop("unknown parameter(s): ", paste(unknown, collapse = ", "), call. = FALSE)
          for (nm in names(pi)) d[[nm]] <- as.numeric(pi[[nm]])
        }
      }
      if (!is.finite(d$N) || d$N < 0 || d$N != floor(d$N) || any(w > d$N)) {
        stop("discrete-binomial N must be an integer >= max(weights)", call. = FALSE)
      }
      cpp_params[[i]] <- c(d$N, d$alpha, d$beta)
      named_params[[i]] <- d
    } else if (tp == "real-normal") {
      if (is_micro) {
        d <- list(m0 = NaN, k0 = NaN, v0 = NaN, nu0 = NaN)
      } else {
        mu <- if (length(w)) mean(w) else 0
        # Match graph-tool's numpy std()^2 (population variance)
        vpop <- if (length(w)) mean((w - mu)^2) else 1
        d <- list(m0 = mu, k0 = 1, v0 = if (vpop > 0) vpop else 1, nu0 = 3)
        if (is.list(pi)) {
          unknown <- setdiff(names(pi), names(d))
          if (length(unknown)) stop("unknown parameter(s): ", paste(unknown, collapse = ", "), call. = FALSE)
          for (nm in names(pi)) d[[nm]] <- as.numeric(pi[[nm]])
        }
      }
      cpp_params[[i]] <- c(d$m0, d$k0, d$v0, d$nu0)
      named_params[[i]] <- d
    }
  }
  list(types = types, cpp_params = cpp_params, params = named_params)
}

new_weighted_state <- function(graph, blocks, weights, spec, degree_corrected,
                               degree_prior, fit = NULL) {
  structure(list(
    blocks = as.integer(blocks),
    graph = graph,
    weights = weights,
    types = spec$types,
    params = spec$params,
    cpp_params = spec$cpp_params,
    degree_corrected = degree_corrected,
    degree_prior = degree_prior,
    num_vertices = graph_num_vertices(graph),
    num_edges = graph_num_edges(graph),
    B = length(unique(blocks)),
    fit = fit
  ), class = "rgraphtool_weighted_block_state")
}

validate_weighted_state <- function(state) {
  if (!inherits(state, "rgraphtool_weighted_block_state")) {
    stop("state must be created by weighted_block_state() or minimize_weighted_blockmodel_dl()",
         call. = FALSE)
  }
  if (graph_num_edges(state$graph) != state$num_edges ||
      graph_num_vertices(state$graph) != state$num_vertices) {
    stop("the graph has changed since the weighted block state was created", call. = FALSE)
  }
}

#' @noRd
#' @export
weighted_block_state <- function(graph, weights, types,
                                 params = "microcanonical", blocks = NULL,
                                 degree_corrected = TRUE,
                                 degree_prior = c("distributed", "uniform")) {
  graph_pointer(graph)
  degree_prior <- match.arg(degree_prior)
  if (!is.logical(degree_corrected) || length(degree_corrected) != 1L ||
      is.na(degree_corrected)) {
    stop("degree_corrected must be TRUE or FALSE", call. = FALSE)
  }
  n <- graph_num_vertices(graph)
  if (is.null(blocks)) {
    blocks <- rep.int(1L, n)
  } else {
    blocks <- validate_label_vector(blocks, "blocks")
    if (length(blocks) != n || anyNA(blocks)) {
      stop("blocks must be an integer vector of length graph_num_vertices(graph)",
           call. = FALSE)
    }
  }
  wmat <- resolve_weight_matrix(graph, weights)
  spec <- resolve_covariate_specs(wmat, types, params)
  new_weighted_state(graph, blocks, wmat, spec, degree_corrected, degree_prior)
}

#' @noRd
#' @export
weighted_description_length <- function(state, rec = TRUE, rec_dl = TRUE,
                                        degree_prior = NULL) {
  validate_weighted_state(state)
  if (is.null(degree_prior)) degree_prior <- state$degree_prior
  degree_prior <- match.arg(degree_prior, c("distributed", "uniform"))
  cpp_weighted_entropy(graph_pointer(state$graph), state$blocks, state$weights,
                       state$types, state$cpp_params, state$degree_corrected,
                       degree_prior == "distributed", isTRUE(rec),
                       isTRUE(rec_dl), q_cache_max())
}

#' @noRd
#' @export
minimize_weighted_blockmodel_dl <- function(graph, weights, types,
                                            params = "microcanonical",
                                            degree_corrected = TRUE,
                                            degree_prior = c("distributed", "uniform"),
                                            B_min = 1L, B_max = Inf,
                                            refine = TRUE, epsilon = 1e-3) {
  pointer <- graph_pointer(graph)
  degree_prior <- match.arg(degree_prior)
  if (!is.logical(degree_corrected) || length(degree_corrected) != 1L ||
      is.na(degree_corrected)) {
    stop("degree_corrected must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.numeric(B_min) || length(B_min) != 1L || is.na(B_min) || B_min < 1 ||
      !is.numeric(B_max) || length(B_max) != 1L || is.na(B_max) || B_max < B_min) {
    stop("B_min and B_max must satisfy 1 <= B_min <= B_max", call. = FALSE)
  }
  wmat <- resolve_weight_matrix(graph, weights)
  spec <- resolve_covariate_specs(wmat, types, params)
  res <- cpp_weighted_minimize(pointer, wmat, spec$types, spec$cpp_params,
                               degree_corrected, degree_prior == "distributed",
                               as.integer(floor(B_min)),
                               if (is.finite(B_max)) floor(B_max) else Inf,
                               isTRUE(refine), as.numeric(epsilon), q_cache_max())
  new_weighted_state(graph, res$blocks, wmat, spec, degree_corrected,
                     degree_prior, fit = list(description_length = res$entropy))
}

#' @noRd
#' @export
print.rgraphtool_weighted_block_state <- function(x, ...) {
  dl <- weighted_description_length(x)
  cat(sprintf("<weighted block state: %s, %s, B = %d, description length %.4f nats>\n",
              if (x$degree_corrected) "degree-corrected" else "non-degree-corrected",
              paste(x$types, collapse = "+"), x$B, dl))
  invisible(x)
}

new_nested_weighted_state <- function(graph, levels, weights, spec,
                                      degree_corrected, degree_prior,
                                      fit = NULL) {
  structure(list(
    levels = levels,
    graph = graph,
    weights = weights,
    types = spec$types,
    params = spec$params,
    cpp_params = spec$cpp_params,
    degree_corrected = degree_corrected,
    degree_prior = degree_prior,
    num_vertices = graph_num_vertices(graph),
    num_edges = graph_num_edges(graph),
    B = vapply(levels, function(b) length(unique(b)), integer(1)),
    fit = fit
  ), class = "rgraphtool_nested_weighted_block_state")
}

validate_nested_weighted_state <- function(state) {
  if (!inherits(state, "rgraphtool_nested_weighted_block_state")) {
    stop("state must be created by nested_weighted_block_state() or minimize_nested_weighted_blockmodel_dl()",
         call. = FALSE)
  }
  if (graph_num_edges(state$graph) != state$num_edges ||
      graph_num_vertices(state$graph) != state$num_vertices) {
    stop("the graph has changed since the nested weighted block state was created",
         call. = FALSE)
  }
}

#' @noRd
#' @export
nested_weighted_block_state <- function(graph, weights, types,
                                        params = "microcanonical",
                                        levels = NULL,
                                        degree_corrected = TRUE,
                                        degree_prior = c("distributed", "uniform")) {
  graph_pointer(graph)
  degree_prior <- match.arg(degree_prior)
  if (!is.logical(degree_corrected) || length(degree_corrected) != 1L ||
      is.na(degree_corrected)) {
    stop("degree_corrected must be TRUE or FALSE", call. = FALSE)
  }
  n <- graph_num_vertices(graph)
  if (n < 1L) stop("nested weighted states need at least one vertex", call. = FALSE)
  if (is.null(levels)) {
    levels <- list(rep.int(1L, n), 1L)
  } else {
    if (!is.list(levels) || length(levels) < 1L) {
      stop("levels must be a non-empty list of 1-based integer label vectors",
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
  wmat <- resolve_weight_matrix(graph, weights)
  spec <- resolve_covariate_specs(wmat, types, params)
  new_nested_weighted_state(graph, levels, wmat, spec, degree_corrected, degree_prior)
}

#' @noRd
#' @export
nested_weighted_level_entropy <- function(state, degree_prior = NULL) {
  validate_nested_weighted_state(state)
  if (is.null(degree_prior)) degree_prior <- state$degree_prior
  degree_prior <- match.arg(degree_prior, c("distributed", "uniform"))
  cpp_nested_weighted_entropy(graph_pointer(state$graph), state$levels,
                              state$weights, state$types, state$cpp_params,
                              state$degree_corrected,
                              degree_prior == "distributed", q_cache_max())
}

#' @noRd
#' @export
nested_weighted_description_length <- function(state, degree_prior = NULL) {
  sum(nested_weighted_level_entropy(state, degree_prior = degree_prior))
}

#' @noRd
#' @export
minimize_nested_weighted_blockmodel_dl <- function(graph, weights, types,
                                                   params = "microcanonical",
                                                   degree_corrected = TRUE,
                                                   degree_prior = c("distributed", "uniform"),
                                                   B_min = 1L, B_max = Inf,
                                                   refine = TRUE, epsilon = 1e-3) {
  pointer <- graph_pointer(graph)
  degree_prior <- match.arg(degree_prior)
  if (graph_num_vertices(graph) < 1L) {
    stop("nested weighted states need at least one vertex", call. = FALSE)
  }
  if (!is.logical(degree_corrected) || length(degree_corrected) != 1L ||
      is.na(degree_corrected)) {
    stop("degree_corrected must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.numeric(B_min) || length(B_min) != 1L || is.na(B_min) || B_min < 1 ||
      !is.numeric(B_max) || length(B_max) != 1L || is.na(B_max) || B_max < B_min) {
    stop("B_min and B_max must satisfy 1 <= B_min <= B_max", call. = FALSE)
  }
  wmat <- resolve_weight_matrix(graph, weights)
  spec <- resolve_covariate_specs(wmat, types, params)
  res <- cpp_nested_weighted_minimize(pointer, wmat, spec$types, spec$cpp_params,
                                      degree_corrected,
                                      degree_prior == "distributed",
                                      as.integer(floor(B_min)),
                                      if (is.finite(B_max)) floor(B_max) else Inf,
                                      isTRUE(refine), as.numeric(epsilon),
                                      q_cache_max())
  new_nested_weighted_state(graph, res$bs, wmat, spec, degree_corrected,
                            degree_prior,
                            fit = list(description_length = res$entropy,
                                       level_entropy = res$level_entropy))
}

#' @noRd
#' @export
print.rgraphtool_nested_weighted_block_state <- function(x, ...) {
  dl <- nested_weighted_description_length(x)
  cat(sprintf("<nested weighted block state: %s, %s, %d levels, B = %s, description length %.4f nats>\n",
              if (x$degree_corrected) "degree-corrected" else "non-degree-corrected",
              paste(x$types, collapse = "+"), length(x$levels),
              paste(x$B, collapse = "-"), dl))
  invisible(x)
}

# Declarative weighted model comparison (with automatic Jacobian adjustment for
# continuous weight transformations such as log-normal vs exponential/normal!)
default_weight_candidates <- function(w) {
  is_int <- all(w >= 0) && all(w == floor(w))
  is_pos <- all(w > 0)
  if (is_int) {
    list(
      "discrete-geometric" = list(type = "discrete-geometric", transform = "identity", params = "microcanonical"),
      "discrete-poisson"   = list(type = "discrete-poisson",   transform = "identity", params = "microcanonical"),
      "discrete-binomial"  = list(type = "discrete-binomial",  transform = "identity", params = "microcanonical")
    )
  } else if (is_pos) {
    list(
      "real-exponential" = list(type = "real-exponential", transform = "identity", params = "microcanonical"),
      "real-normal"      = list(type = "real-normal",      transform = "identity", params = "microcanonical"),
      "log-normal"       = list(type = "real-normal",      transform = "log",      params = "microcanonical")
    )
  } else {
    list(
      "real-normal (micro)" = list(type = "real-normal", transform = "identity", params = "microcanonical"),
      "real-normal (bayes)" = list(type = "real-normal", transform = "identity", params = "bayesian")
    )
  }
}

#' @noRd
#' @export
compare_weight_models <- function(graph, weights, candidates = NULL,
                                  nested = TRUE, degree_corrected = TRUE,
                                  degree_prior = c("distributed", "uniform"),
                                  method = c("map", "posterior"),
                                  evidence_method = c("variational", "importance"),
                                  nruns = 1L, nsamples = 200L, thin = 5L,
                                  burnin = NULL, mle = TRUE,
                                  B_min = 1L, B_max = Inf, refine = TRUE) {
  graph_pointer(graph)
  degree_prior <- match.arg(degree_prior)
  method <- match.arg(method)
  evidence_method <- match.arg(evidence_method)
  if (!is.logical(nested) || length(nested) != 1L || is.na(nested)) {
    stop("nested must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.numeric(nruns) || length(nruns) != 1L || is.na(nruns) ||
      nruns < 1 || nruns != floor(nruns)) {
    stop("nruns must be a positive integer", call. = FALSE)
  }
  nruns <- as.integer(nruns)

  raw_mat <- resolve_weight_matrix(graph, weights)
  if (ncol(raw_mat) != 1L) {
    stop("compare_weight_models() expects a single edge weight covariate", call. = FALSE)
  }
  w_raw <- raw_mat[, 1L]

  if (is.null(candidates)) {
    candidates <- default_weight_candidates(w_raw)
  } else if (is.character(candidates)) {
    # Shorthand character names including "log-normal"
    cands <- list()
    for (nm in candidates) {
      if (nm == "log-normal") {
        cands[[nm]] <- list(type = "real-normal", transform = "log", params = "microcanonical")
      } else if (nm %in% valid_covariate_types()) {
        cands[[nm]] <- list(type = nm, transform = "identity", params = "microcanonical")
      } else {
        stop("unknown weight candidate: '", nm, "'", call. = FALSE)
      }
    }
    candidates <- cands
  } else if (!is.list(candidates) || is.null(names(candidates)) ||
             any(!nzchar(names(candidates)))) {
    stop("candidates must be NULL, a character vector, or a named list of candidate specs",
         call. = FALSE)
  }

  names_c <- names(candidates)
  states <- stats::setNames(vector("list", length(names_c)), names_c)
  evidences <- if (method == "posterior") {
    stats::setNames(vector("list", length(names_c)), names_c)
  } else {
    NULL
  }
  rows <- vector("list", length(names_c))

  for (i in seq_along(names_c)) {
    nm <- names_c[i]
    sp <- candidates[[nm]]
    tp <- sp$type
    tr <- if (is.null(sp$transform)) "identity" else sp$transform
    sh <- if (is.null(sp$shift)) 0 else sp$shift
    pr <- if (is.null(sp$params)) "microcanonical" else sp$params

    w_tr <- transform_weights(w_raw, transform = tr, shift = sh)
    log_jac <- attr(w_tr, "log_jacobian")
    # Change of variables: p_W(w) = p_Y(f(w)) * |f'(w)| =>
    # -log p_W(w) = -log p_Y(f(w)) - sum log |f'(w)| = DL_Y - log_jac
    jac_adj <- if (is.finite(log_jac)) -log_jac else 0

    best_st <- NULL
    best_dl <- Inf
    for (run in seq_len(nruns)) {
      st <- if (nested) {
        minimize_nested_weighted_blockmodel_dl(
          graph, weights = as.numeric(w_tr), types = tp, params = pr,
          degree_corrected = degree_corrected, degree_prior = degree_prior,
          B_min = B_min, B_max = B_max, refine = refine
        )
      } else {
        minimize_weighted_blockmodel_dl(
          graph, weights = as.numeric(w_tr), types = tp, params = pr,
          degree_corrected = degree_corrected, degree_prior = degree_prior,
          B_min = B_min, B_max = B_max, refine = refine
        )
      }
      dl_raw <- if (nested) nested_weighted_description_length(st) else weighted_description_length(st)
      if (dl_raw < best_dl) {
        best_dl <- dl_raw
        best_st <- st
      }
    }
    states[[nm]] <- best_st
    dl_adj <- best_dl + jac_adj

    if (method == "posterior") {
      ev <- model_evidence(best_st, nsamples = nsamples, thin = thin,
                           burnin = burnin, mle = mle,
                           method = evidence_method,
                           B_min = B_min, B_max = B_max)
      evidences[[nm]] <- ev
      neg_log_ev <- ev$neg_log_evidence + jac_adj
      post_h <- ev$posterior_entropy
    } else {
      neg_log_ev <- dl_adj
      post_h <- NA_real_
    }

    rows[[i]] <- data.frame(
      model = nm,
      type = tp,
      transform = tr,
      B = if (nested) best_st$B[1L] else best_st$B,
      raw_dl = best_dl,
      jacobian_adj = jac_adj,
      description_length = dl_adj,
      neg_log_evidence = neg_log_ev,
      posterior_entropy = post_h,
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
    for (j in seq_len(k)) {
      lbf[i, j] <- scores[j] - scores[i]
    }
  }

  best_model <- tab$model[1L]
  structure(list(
    summary = tab,
    best_model = best_model,
    best_state = states[[best_model]],
    states = states[tab$model],
    evidences = if (is.null(evidences)) NULL else evidences[tab$model],
    log_bayes_factors = lbf,
    nested = nested,
    degree_corrected = degree_corrected,
    degree_prior = degree_prior,
    method = method,
    evidence_method = if (method == "posterior") evidence_method else NULL
  ), class = "rgraphtool_weight_model_comparison")
}

#' @noRd
#' @export
print.rgraphtool_weight_model_comparison <- function(x, ...) {
  cat(sprintf("<weighted SBM model comparison (%s, %s, %s)>\n",
              if (x$nested) "nested" else "flat",
              if (x$degree_corrected) "DC" else "NDC",
              if (x$method == "posterior") paste0("posterior/", x$evidence_method) else "MAP"))
  df <- x$summary
  cols <- if (x$method == "posterior") {
    c("model", "B", "description_length", "jacobian_adj", "posterior_entropy",
      "neg_log_evidence", "delta", "posterior_weight")
  } else {
    c("model", "B", "description_length", "jacobian_adj", "delta", "posterior_weight")
  }
  out <- df[, cols, drop = FALSE]
  for (col in intersect(c("description_length", "jacobian_adj",
                          "posterior_entropy", "neg_log_evidence", "delta"),
                        names(out))) {
    out[[col]] <- round(out[[col]], 3)
  }
  out$posterior_weight <- signif(out$posterior_weight, 3)
  print(out, row.names = FALSE)
  cat(sprintf("Best weight model: %s\n", x$best_model))
  invisible(x)
}
