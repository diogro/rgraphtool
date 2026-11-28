resolve_layer_vector <- function(graph, layers) {
  e <- graph_num_edges(graph)
  raw <- NULL
  if (is.character(layers) && length(layers) == 1L && e != 1L) {
    if (is.na(layers) || !nzchar(layers)) {
      stop("layers property name must be a non-empty character string", call. = FALSE)
    }
    raw <- graph_edge_property(graph, layers)
  } else if (is.character(layers) && length(layers) == 1L && e == 1L &&
             layers %in% graph_edge_properties(graph)) {
    raw <- graph_edge_property(graph, layers)
  } else if (is.numeric(layers) || is.integer(layers) ||
             is.character(layers) || is.factor(layers) || is.logical(layers)) {
    if (length(layers) != e) {
      stop("layers vector must have length graph_num_edges(graph)", call. = FALSE)
    }
    raw <- layers
  } else {
    stop("layers must be an edge property name or a vector of length graph_num_edges(graph)",
         call. = FALSE)
  }
  if (anyNA(raw)) {
    stop("layers must not contain missing values", call. = FALSE)
  }
  if (is.numeric(raw) && any(!is.finite(raw))) {
    stop("numeric layers must be finite", call. = FALSE)
  }
  if (e == 0L) {
    return(list(codes = integer(0), levels = character(0), C = 0L))
  }
  f <- factor(raw)
  codes <- as.integer(f)
  lvls <- levels(f)
  list(codes = codes, levels = lvls, C = length(lvls))
}

resolve_layered_weights <- function(graph, weights = NULL, types = NULL,
                                    params = "microcanonical") {
  if (is.null(weights)) {
    return(list(
      W = matrix(0, nrow = 0L, ncol = 0L),
      types = character(0),
      params = list(),
      cpp_params = list(),
      weighted = FALSE
    ))
  }
  if (is.null(types)) {
    stop("types must be specified when weights is provided", call. = FALSE)
  }
  W <- resolve_weight_matrix(graph, weights)
  spec <- resolve_covariate_specs(W, types, params = params)
  list(
    W = W,
    types = spec$types,
    params = spec$params,
    cpp_params = spec$cpp_params,
    weighted = TRUE
  )
}

new_layered_state <- function(graph, blocks, layer_info, independent,
                              degree_corrected, degree_prior, fit = NULL,
                              wspec = NULL) {
  if (is.null(wspec)) {
    wspec <- list(
      W = matrix(0, nrow = 0L, ncol = 0L),
      types = character(0),
      params = list(),
      cpp_params = list(),
      weighted = FALSE
    )
  }
  structure(list(
    blocks = as.integer(blocks),
    graph = graph,
    layers = layer_info$codes,
    layer_levels = layer_info$levels,
    C = layer_info$C,
    independent = isTRUE(independent),
    degree_corrected = isTRUE(degree_corrected),
    degree_prior = degree_prior,
    weights = wspec$W,
    types = wspec$types,
    params = wspec$params,
    cpp_params = wspec$cpp_params,
    weighted = isTRUE(wspec$weighted),
    num_vertices = graph_num_vertices(graph),
    num_edges = graph_num_edges(graph),
    B = length(unique(blocks)),
    fit = fit
  ), class = "rgraphtool_layered_block_state")
}

validate_layered_state <- function(state) {
  if (!inherits(state, "rgraphtool_layered_block_state")) {
    stop("state must be created by layered_block_state() or minimize_layered_blockmodel_dl()",
         call. = FALSE)
  }
  if (graph_num_edges(state$graph) != state$num_edges ||
      graph_num_vertices(state$graph) != state$num_vertices) {
    stop("the graph has changed since the layered block state was created", call. = FALSE)
  }
}

#' @noRd
#' @export
layered_block_state <- function(graph, layers, independent = TRUE,
                                blocks = NULL, degree_corrected = TRUE,
                                degree_prior = c("distributed", "uniform"),
                                weights = NULL, types = NULL,
                                params = "microcanonical") {
  graph_pointer(graph)
  degree_prior <- match.arg(degree_prior)
  if (!is.logical(independent) || length(independent) != 1L || is.na(independent)) {
    stop("independent must be TRUE or FALSE", call. = FALSE)
  }
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
  linfo <- resolve_layer_vector(graph, layers)
  wspec <- resolve_layered_weights(graph, weights, types, params)
  new_layered_state(graph, blocks, linfo, independent, degree_corrected,
                    degree_prior, wspec = wspec)
}

#' @noRd
#' @export
layered_description_length <- function(state, lmemb = TRUE, degree_prior = NULL) {
  validate_layered_state(state)
  if (is.null(degree_prior)) degree_prior <- state$degree_prior
  degree_prior <- match.arg(degree_prior, c("distributed", "uniform"))
  w_mat <- if (isTRUE(state$weighted)) state$weights else matrix(0, nrow = 0L, ncol = 0L)
  w_types <- if (isTRUE(state$weighted)) state$types else character(0)
  w_params <- if (isTRUE(state$weighted)) state$cpp_params else list()
  cpp_layered_entropy(graph_pointer(state$graph), state$blocks, state$layers,
                      state$independent, state$degree_corrected,
                      degree_prior == "distributed", isTRUE(lmemb), q_cache_max(),
                      w_mat, w_types, w_params)
}

#' @noRd
#' @export
minimize_layered_blockmodel_dl <- function(graph, layers, independent = TRUE,
                                           degree_corrected = TRUE,
                                           degree_prior = c("distributed", "uniform"),
                                           B_min = 1L, B_max = Inf,
                                           refine = TRUE, epsilon = 1e-3,
                                           weights = NULL, types = NULL,
                                           params = "microcanonical",
                                           nthreads = NULL) {
  pointer <- graph_pointer(graph)
  degree_prior <- match.arg(degree_prior)
  if (!is.logical(independent) || length(independent) != 1L || is.na(independent)) {
    stop("independent must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(degree_corrected) || length(degree_corrected) != 1L ||
      is.na(degree_corrected)) {
    stop("degree_corrected must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.numeric(B_min) || length(B_min) != 1L || is.na(B_min) || B_min < 1 ||
      !is.numeric(B_max) || length(B_max) != 1L || is.na(B_max) || B_max < B_min) {
    stop("B_min and B_max must satisfy 1 <= B_min <= B_max", call. = FALSE)
  }
  nt <- if (is.null(nthreads)) 0L else as.integer(nthreads)
  linfo <- resolve_layer_vector(graph, layers)
  wspec <- resolve_layered_weights(graph, weights, types, params)
  res <- cpp_layered_minimize(pointer, linfo$codes, independent,
                              degree_corrected, degree_prior == "distributed",
                              as.integer(floor(B_min)),
                              if (is.finite(B_max)) floor(B_max) else Inf,
                              isTRUE(refine), as.numeric(epsilon), q_cache_max(),
                              wspec$W, wspec$types, wspec$cpp_params, nt)
  new_layered_state(graph, res$blocks, linfo, independent, degree_corrected,
                    degree_prior, fit = list(description_length = res$entropy),
                    wspec = wspec)
}

#' @noRd
#' @export
print.rgraphtool_layered_block_state <- function(x, ...) {
  dl <- layered_description_length(x)
  wtxt <- if (isTRUE(x$weighted)) paste0(", ", paste(x$types, collapse = "+")) else ""
  cat(sprintf("<layered block state: %s, %d %s layer%s%s, B = %d, description length %.4f nats>\n",
              if (x$degree_corrected) "degree-corrected" else "non-degree-corrected",
              x$C,
              if (x$independent) "independent" else "edge-covariate",
              if (x$C == 1L) "" else "s",
              wtxt, x$B, dl))
  invisible(x)
}

new_nested_layered_state <- function(graph, levels, layer_info, independent,
                                     degree_corrected, degree_prior,
                                     fit = NULL, wspec = NULL) {
  if (is.null(wspec)) {
    wspec <- list(
      W = matrix(0, nrow = 0L, ncol = 0L),
      types = character(0),
      params = list(),
      cpp_params = list(),
      weighted = FALSE
    )
  }
  structure(list(
    levels = levels,
    graph = graph,
    layers = layer_info$codes,
    layer_levels = layer_info$levels,
    C = layer_info$C,
    independent = isTRUE(independent),
    degree_corrected = isTRUE(degree_corrected),
    degree_prior = degree_prior,
    weights = wspec$W,
    types = wspec$types,
    params = wspec$params,
    cpp_params = wspec$cpp_params,
    weighted = isTRUE(wspec$weighted),
    num_vertices = graph_num_vertices(graph),
    num_edges = graph_num_edges(graph),
    B = vapply(levels, function(b) length(unique(b)), integer(1)),
    fit = fit
  ), class = "rgraphtool_nested_layered_block_state")
}

validate_nested_layered_state <- function(state) {
  if (!inherits(state, "rgraphtool_nested_layered_block_state")) {
    stop("state must be created by nested_layered_block_state() or minimize_nested_layered_blockmodel_dl()",
         call. = FALSE)
  }
  if (graph_num_edges(state$graph) != state$num_edges ||
      graph_num_vertices(state$graph) != state$num_vertices) {
    stop("the graph has changed since the nested layered block state was created",
         call. = FALSE)
  }
}

#' @noRd
#' @export
nested_layered_block_state <- function(graph, layers, independent = TRUE,
                                       levels = NULL, degree_corrected = TRUE,
                                       degree_prior = c("distributed", "uniform"),
                                       weights = NULL, types = NULL,
                                       params = "microcanonical") {
  graph_pointer(graph)
  degree_prior <- match.arg(degree_prior)
  if (!is.logical(independent) || length(independent) != 1L || is.na(independent)) {
    stop("independent must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(degree_corrected) || length(degree_corrected) != 1L ||
      is.na(degree_corrected)) {
    stop("degree_corrected must be TRUE or FALSE", call. = FALSE)
  }
  n <- graph_num_vertices(graph)
  if (n < 1L) stop("nested layered states need at least one vertex", call. = FALSE)
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
  linfo <- resolve_layer_vector(graph, layers)
  wspec <- resolve_layered_weights(graph, weights, types, params)
  new_nested_layered_state(graph, levels, linfo, independent, degree_corrected,
                           degree_prior, wspec = wspec)
}

#' @noRd
#' @export
nested_layered_level_entropy <- function(state, degree_prior = NULL) {
  validate_nested_layered_state(state)
  if (is.null(degree_prior)) degree_prior <- state$degree_prior
  degree_prior <- match.arg(degree_prior, c("distributed", "uniform"))
  w_mat <- if (isTRUE(state$weighted)) state$weights else matrix(0, nrow = 0L, ncol = 0L)
  w_types <- if (isTRUE(state$weighted)) state$types else character(0)
  w_params <- if (isTRUE(state$weighted)) state$cpp_params else list()
  cpp_nested_layered_entropy(graph_pointer(state$graph), state$levels,
                             state$layers, state$independent,
                             state$degree_corrected,
                             degree_prior == "distributed", q_cache_max(),
                             w_mat, w_types, w_params)
}

#' @noRd
#' @export
nested_layered_description_length <- function(state, degree_prior = NULL) {
  sum(nested_layered_level_entropy(state, degree_prior = degree_prior))
}

#' @noRd
#' @export
minimize_nested_layered_blockmodel_dl <- function(graph, layers, independent = TRUE,
                                                  degree_corrected = TRUE,
                                                  degree_prior = c("distributed", "uniform"),
                                                  B_min = 1L, B_max = Inf,
                                                  refine = TRUE, epsilon = 1e-3,
                                                  weights = NULL, types = NULL,
                                                  params = "microcanonical",
                                                  nthreads = NULL) {
  pointer <- graph_pointer(graph)
  degree_prior <- match.arg(degree_prior)
  if (graph_num_vertices(graph) < 1L) {
    stop("nested layered states need at least one vertex", call. = FALSE)
  }
  if (!is.logical(independent) || length(independent) != 1L || is.na(independent)) {
    stop("independent must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(degree_corrected) || length(degree_corrected) != 1L ||
      is.na(degree_corrected)) {
    stop("degree_corrected must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.numeric(B_min) || length(B_min) != 1L || is.na(B_min) || B_min < 1 ||
      !is.numeric(B_max) || length(B_max) != 1L || is.na(B_max) || B_max < B_min) {
    stop("B_min and B_max must satisfy 1 <= B_min <= B_max", call. = FALSE)
  }
  nt <- if (is.null(nthreads)) 0L else as.integer(nthreads)
  linfo <- resolve_layer_vector(graph, layers)
  wspec <- resolve_layered_weights(graph, weights, types, params)
  res <- cpp_nested_layered_minimize(pointer, linfo$codes, independent,
                                     degree_corrected,
                                     degree_prior == "distributed",
                                     as.integer(floor(B_min)),
                                     if (is.finite(B_max)) floor(B_max) else Inf,
                                     isTRUE(refine), as.numeric(epsilon),
                                     q_cache_max(), wspec$W, wspec$types,
                                     wspec$cpp_params, nt)
  new_nested_layered_state(graph, res$bs, linfo, independent, degree_corrected,
                           degree_prior,
                           fit = list(description_length = res$entropy,
                                      level_entropy = res$level_entropy),
                           wspec = wspec)
}

#' @noRd
#' @export
print.rgraphtool_nested_layered_block_state <- function(x, ...) {
  dl <- nested_layered_description_length(x)
  wtxt <- if (isTRUE(x$weighted)) paste0(", ", paste(x$types, collapse = "+")) else ""
  cat(sprintf("<nested layered block state: %s, %d %s layer%s%s, %d levels, B = %s, description length %.4f nats>\n",
              if (x$degree_corrected) "degree-corrected" else "non-degree-corrected",
              x$C,
              if (x$independent) "independent" else "edge-covariate",
              if (x$C == 1L) "" else "s",
              wtxt, length(x$levels), paste(x$B, collapse = "-"), dl))
  invisible(x)
}

# Declarative layered model comparison (edge-covariate vs independent-layers,
# optionally across DC/NDC and MAP vs full Bayesian posterior marginal evidence)
#' @noRd
#' @export
compare_layered_models <- function(graph, layers,
                                   models = c("edge_covariates", "independent_layers"),
                                   nested = TRUE, degree_corrected = TRUE,
                                   degree_prior = c("distributed", "uniform"),
                                   method = c("map", "posterior"),
                                   evidence_method = c("variational", "importance"),
                                   nruns = 1L, cores = NULL, nsamples = 200L, thin = 5L,
                                   burnin = NULL, mle = TRUE,
                                   B_min = 1L, B_max = Inf, refine = TRUE,
                                   weights = NULL, types = NULL,
                                   params = "microcanonical",
                                   nthreads = NULL) {
  graph_pointer(graph)
  degree_prior <- match.arg(degree_prior)
  method <- match.arg(method)
  evidence_method <- match.arg(evidence_method)
  if (!is.logical(nested) || length(nested) != 1L || is.na(nested)) {
    stop("nested must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(degree_corrected) || length(degree_corrected) != 1L ||
      is.na(degree_corrected)) {
    stop("degree_corrected must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.numeric(nruns) || length(nruns) != 1L || is.na(nruns) ||
      nruns < 1 || nruns != floor(nruns)) {
    stop("nruns must be a positive integer", call. = FALSE)
  }
  nruns <- as.integer(nruns)

  valid_models <- c("edge_covariates", "independent_layers",
                    "edge_covariates_dc", "edge_covariates_ndc",
                    "independent_layers_dc", "independent_layers_ndc",
                    "pp_edge_covariates", "pp_independent_layers",
                    "pp_edge_covariates_uniform",
                    "pp_independent_layers_uniform")
  if (!is.character(models) || !length(models) || anyNA(models) ||
      any(!models %in% valid_models)) {
    stop("models must be a non-empty subset of: ",
         paste(valid_models, collapse = ", "), call. = FALSE)
  }
  models <- unique(models)

  # Resolve layers & optional weights once for validation
  linfo <- resolve_layer_vector(graph, layers)
  wspec <- resolve_layered_weights(graph, weights, types, params)

  spec_for <- function(nm) {
    switch(
      nm,
      edge_covariates        = list(independent = FALSE, dc = degree_corrected, pp = FALSE, uniform = FALSE),
      independent_layers     = list(independent = TRUE,  dc = degree_corrected, pp = FALSE, uniform = FALSE),
      edge_covariates_dc     = list(independent = FALSE, dc = TRUE,  pp = FALSE, uniform = FALSE),
      edge_covariates_ndc    = list(independent = FALSE, dc = FALSE, pp = FALSE, uniform = FALSE),
      independent_layers_dc  = list(independent = TRUE,  dc = TRUE,  pp = FALSE, uniform = FALSE),
      independent_layers_ndc = list(independent = TRUE,  dc = FALSE, pp = FALSE, uniform = FALSE),
      pp_edge_covariates         = list(independent = FALSE, dc = degree_corrected, pp = TRUE, uniform = FALSE),
      pp_independent_layers      = list(independent = TRUE,  dc = degree_corrected, pp = TRUE, uniform = FALSE),
      pp_edge_covariates_uniform = list(independent = FALSE, dc = degree_corrected, pp = TRUE, uniform = TRUE),
      pp_independent_layers_uniform = list(independent = TRUE, dc = degree_corrected, pp = TRUE, uniform = TRUE)
    )
  }

  states <- stats::setNames(vector("list", length(models)), models)
  evidences <- if (method == "posterior") {
    stats::setNames(vector("list", length(models)), models)
  } else {
    NULL
  }
  rows <- vector("list", length(models))

  for (i in seq_along(models)) {
    nm <- models[i]
    sp <- spec_for(nm)
    best_st <- NULL
    best_dl <- Inf
    .starts <- parallel_starts(nruns, function(run) {
      st <- if (isTRUE(sp$pp)) {
        minimize_pp_blockmodel_dl(
          graph, uniform = sp$uniform, degree_corrected = sp$dc,
          degree_prior = degree_prior, B_min = B_min, B_max = B_max,
          refine = refine,
          weights = if (wspec$weighted) wspec$W else NULL,
          types = if (wspec$weighted) wspec$types else NULL,
          params = if (wspec$weighted) params else "microcanonical",
          layers = linfo$codes, independent = sp$independent
        )
      } else if (nested) {
        minimize_nested_layered_blockmodel_dl(
          graph, layers = linfo$codes, independent = sp$independent,
          degree_corrected = sp$dc, degree_prior = degree_prior,
          B_min = B_min, B_max = B_max, refine = refine,
          weights = if (wspec$weighted) wspec$W else NULL,
          types = if (wspec$weighted) wspec$types else NULL,
          params = if (wspec$weighted) params else "microcanonical",
          nthreads = nthreads
        )
      } else {
        minimize_layered_blockmodel_dl(
          graph, layers = linfo$codes, independent = sp$independent,
          degree_corrected = sp$dc, degree_prior = degree_prior,
          B_min = B_min, B_max = B_max, refine = refine,
          weights = if (wspec$weighted) wspec$W else NULL,
          types = if (wspec$weighted) wspec$types else NULL,
          params = if (wspec$weighted) params else "microcanonical",
          nthreads = nthreads
        )
      }
      st
    }, cores = cores, graph = graph)
    for (run in seq_len(nruns)) {
      st <- .starts[[run]]
      st$layer_levels <- linfo$levels
      dl <- state_description_length(st)
      if (dl < best_dl) {
        best_dl <- dl
        best_st <- st
      }
    }
    states[[nm]] <- best_st

    if (method == "posterior") {
      ev <- model_evidence(best_st, nsamples = nsamples, thin = thin,
                           burnin = burnin, mle = mle,
                           method = evidence_method,
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
      independent = sp$independent,
      degree_corrected = sp$dc,
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
    C = linfo$C,
    layer_levels = linfo$levels,
    degree_prior = degree_prior,
    method = method,
    evidence_method = if (method == "posterior") evidence_method else NULL
  ), class = "rgraphtool_layered_model_comparison")
}

#' @noRd
#' @export
print.rgraphtool_layered_model_comparison <- function(x, ...) {
  cat(sprintf("<layered SBM model comparison (%s, %d layers, %s)>\n",
              if (x$nested) "nested" else "flat",
              x$C,
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
  cat(sprintf("Best layered model: %s\n", x$best_model))
  invisible(x)
}
