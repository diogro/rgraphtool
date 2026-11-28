new_pp_state <- function(graph, blocks, uniform, degree_corrected, degree_prior,
                         fit = NULL, wspec = NULL, linfo = NULL,
                         independent = FALSE, present = NULL) {
  if (is.null(wspec)) {
    wspec <- list(
      W = matrix(0, nrow = 0L, ncol = 0L),
      types = character(0),
      params = list(),
      cpp_params = list(),
      weighted = FALSE
    )
  }
  layered <- !is.null(linfo)
  structure(list(
    blocks = as.integer(blocks),
    graph = graph,
    uniform = isTRUE(uniform),
    degree_corrected = isTRUE(degree_corrected),
    degree_prior = degree_prior,
    layers = if (layered) linfo$codes else integer(0),
    layer_levels = if (layered) linfo$levels else character(0),
    C = if (layered) linfo$C else 1L,
    independent = if (layered) isTRUE(independent) else FALSE,
    layered = layered,
    weights = wspec$W,
    types = wspec$types,
    params = wspec$params,
    cpp_params = wspec$cpp_params,
    weighted = isTRUE(wspec$weighted),
    present = present,
    num_vertices = graph_num_vertices(graph),
    num_edges = graph_num_edges(graph),
    B = length(unique(blocks)),
    fit = fit
  ), class = "rgraphtool_pp_block_state")
}

validate_pp_state <- function(state) {
  if (!inherits(state, "rgraphtool_pp_block_state")) {
    stop("state must be created by pp_block_state() or minimize_pp_blockmodel_dl()",
         call. = FALSE)
  }
  if (graph_num_edges(state$graph) != state$num_edges ||
      graph_num_vertices(state$graph) != state$num_vertices) {
    stop("the graph has changed since the planted-partition state was created",
         call. = FALSE)
  }
}

#' @noRd
#' @export
pp_block_state <- function(graph, blocks = NULL, uniform = FALSE,
                           degree_corrected = TRUE,
                           degree_prior = c("distributed", "uniform"),
                           weights = NULL, types = NULL,
                           params = "microcanonical",
                           layers = NULL, independent = FALSE, present = NULL) {
  graph_pointer(graph)
  degree_prior <- match.arg(degree_prior)
  if (!is.logical(uniform) || length(uniform) != 1L || is.na(uniform)) {
    stop("uniform must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(degree_corrected) || length(degree_corrected) != 1L ||
      is.na(degree_corrected)) {
    stop("degree_corrected must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(independent) || length(independent) != 1L || is.na(independent)) {
    stop("independent must be TRUE or FALSE", call. = FALSE)
  }
  wspec <- resolve_layered_weights(graph, weights, types, params)
  linfo <- if (!is.null(layers)) resolve_layer_vector(graph, layers) else NULL
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
  pres <- if (!is.null(present)) {
    if (is.null(linfo)) stop("present needs layers", call. = FALSE)
    resolve_presence(graph, present, linfo)
  } else {
    NULL
  }
  new_pp_state(graph, blocks, uniform, degree_corrected, degree_prior,
               wspec = wspec, linfo = linfo, independent = independent,
               present = pres)
}

#' @noRd
#' @export
pp_description_length <- function(state, degree_prior = NULL) {
  validate_pp_state(state)
  if (is.null(degree_prior)) degree_prior <- state$degree_prior
  degree_prior <- match.arg(degree_prior, c("distributed", "uniform"))
  w_mat <- if (isTRUE(state$weighted)) state$weights else NULL
  w_types <- if (isTRUE(state$weighted)) state$types else NULL
  w_params <- if (isTRUE(state$weighted)) state$cpp_params else NULL
  l_vec <- if (isTRUE(state$layered)) state$layers else NULL
  cpp_pp_entropy(graph_pointer(state$graph), state$blocks, state$uniform,
                 state$degree_corrected, degree_prior == "distributed",
                 q_cache_max(), w_mat, w_types, w_params,
                 l_vec, isTRUE(state$independent), state$present)
}

#' @noRd
#' @export
pp_assortativity <- function(state) {
  validate_pp_state(state)
  w_mat <- if (isTRUE(state$weighted)) state$weights else NULL
  w_types <- if (isTRUE(state$weighted)) state$types else NULL
  w_params <- if (isTRUE(state$weighted)) state$cpp_params else NULL
  l_vec <- if (isTRUE(state$layered)) state$layers else NULL
  st <- cpp_pp_stats(graph_pointer(state$graph), state$blocks, state$uniform,
                     state$degree_corrected,
                     identical(state$degree_prior, "distributed"),
                     q_cache_max(), w_mat, w_types, w_params,
                     l_vec, isTRUE(state$independent), state$present)
  sz <- as.numeric(tabulate(state$blocks, nbins = max(state$blocks, 0L)))
  sz <- sz[sz > 0]
  p <- sz / sum(sz)
  Be <- if (length(p)) exp(-sum(p * log(p))) else 0
  frac <- if (st$E_total > 0) st$E_in / st$E_total else NA_real_
  out <- list(
    E_in = st$E_in,
    E_out = st$E_out,
    E_total = st$E_total,
    fraction_in = frac,
    B = as.integer(st$B),
    Be = Be,
    err = st$err
  )
  if (isTRUE(state$layered)) {
    lnms <- if (length(state$layer_levels) == length(st$E_total_layer)) {
      state$layer_levels
    } else {
      paste0("layer", seq_along(st$E_total_layer))
    }
    out$E_in_layer <- stats::setNames(as.numeric(st$E_in_layer), lnms)
    out$E_out_layer <- stats::setNames(as.numeric(st$E_out_layer), lnms)
    out$E_total_layer <- stats::setNames(as.numeric(st$E_total_layer), lnms)
    out$fraction_in_layer <- ifelse(out$E_total_layer > 0,
                                    out$E_in_layer / out$E_total_layer, NA_real_)
    err_mat <- st$err_layer
    colnames(err_mat) <- lnms
    out$err_layer <- err_mat
  }
  if (isTRUE(state$weighted)) {
    out$weight_in_mean <- as.numeric(st$weight_in_mean)
    out$weight_out_mean <- as.numeric(st$weight_out_mean)
    if (isTRUE(state$layered)) {
      lnms <- names(out$E_total_layer)
      win_m <- st$weight_in_mean_layer
      wout_m <- st$weight_out_mean_layer
      rownames(win_m) <- lnms
      rownames(wout_m) <- lnms
      out$weight_in_mean_layer <- win_m
      out$weight_out_mean_layer <- wout_m
    }
  }
  out
}

#' @noRd
#' @export
minimize_pp_blockmodel_dl <- function(graph, uniform = FALSE,
                                      degree_corrected = TRUE,
                                      degree_prior = c("distributed", "uniform"),
                                      B_min = 1L, B_max = Inf,
                                      refine = TRUE, epsilon = 1e-3,
                                      weights = NULL, types = NULL,
                                      params = "microcanonical",
                                      layers = NULL, independent = FALSE,
                                      present = NULL, spectral_starts = TRUE) {
  pointer <- graph_pointer(graph)
  degree_prior <- match.arg(degree_prior)
  if (!is.logical(uniform) || length(uniform) != 1L || is.na(uniform)) {
    stop("uniform must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(degree_corrected) || length(degree_corrected) != 1L ||
      is.na(degree_corrected)) {
    stop("degree_corrected must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(independent) || length(independent) != 1L || is.na(independent)) {
    stop("independent must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.numeric(B_min) || length(B_min) != 1L || is.na(B_min) || B_min < 1 ||
      !is.numeric(B_max) || length(B_max) != 1L || is.na(B_max) || B_max < B_min) {
    stop("B_min and B_max must satisfy 1 <= B_min <= B_max", call. = FALSE)
  }
  wspec <- resolve_layered_weights(graph, weights, types, params)
  linfo <- if (!is.null(layers)) resolve_layer_vector(graph, layers) else NULL
  w_mat <- if (wspec$weighted) wspec$W else NULL
  w_types <- if (wspec$weighted) wspec$types else NULL
  w_params <- if (wspec$weighted) wspec$cpp_params else NULL
  l_vec <- if (!is.null(linfo)) linfo$codes else NULL
  pres <- if (!is.null(present) && !is.null(linfo)) resolve_presence(graph, present, linfo) else NULL
  res <- cpp_pp_minimize(pointer, uniform, degree_corrected,
                         degree_prior == "distributed",
                         as.integer(floor(B_min)),
                         if (is.finite(B_max)) floor(B_max) else Inf,
                         isTRUE(refine), as.numeric(epsilon), q_cache_max(),
                         w_mat, w_types, w_params, l_vec, isTRUE(independent), pres,
                         if (isTRUE(spectral_starts) && wspec$weighted)
                           weight_spectral_seeds(graph, w_mat, w_types) else NULL)
  new_pp_state(graph, res$blocks, uniform, degree_corrected, degree_prior,
               fit = list(description_length = res$entropy),
               wspec = wspec, linfo = linfo, independent = independent,
               present = pres)
}

#' @noRd
#' @export
print.rgraphtool_pp_block_state <- function(x, ...) {
  dl <- pp_description_length(x)
  asrt <- pp_assortativity(x)
  ltxt <- if (isTRUE(x$layered)) {
    sprintf(", %d layers (%s)", x$C,
            if (isTRUE(x$independent)) "independent" else "edge-covariates")
  } else {
    ""
  }
  wtxt <- if (isTRUE(x$weighted)) {
    sprintf(", weighted (%s)", paste(x$types, collapse = ", "))
  } else {
    ""
  }
  cat(sprintf("<planted-partition block state: %s, %s assortativity%s%s, B = %d (Be = %.2f), E_in/E = %.1f%%, DL = %.4f nats>\n",
              if (x$degree_corrected) "degree-corrected" else "non-degree-corrected",
              if (x$uniform) "uniform" else "non-uniform",
              ltxt,
              wtxt,
              x$B, asrt$Be,
              if (is.finite(asrt$fraction_in)) 100 * asrt$fraction_in else 0,
              dl))
  invisible(x)
}

# Declarative comparison of assortative planted-partition vs general SBM
#' @noRd
#' @export
compare_assortative_models <- function(graph,
                                       models = c("pp", "pp_uniform",
                                                  "sbm_flat", "sbm_nested"),
                                       degree_corrected = TRUE,
                                       degree_prior = c("distributed", "uniform"),
                                       method = c("map", "posterior"),
                                       evidence_method = c("variational", "importance"),
                                       nruns = 1L, cores = NULL, nsamples = 200L, thin = 5L,
                                       burnin = NULL, mle = TRUE,
                                       B_min = 1L, B_max = Inf, refine = TRUE,
                                       weights = NULL, types = NULL,
                                       params = "microcanonical") {
  graph_pointer(graph)
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
  wspec <- resolve_layered_weights(graph, weights, types, params)
  valid <- c("pp", "pp_uniform", "sbm_flat", "sbm_nested")
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
    .starts <- parallel_starts(nruns, function(run) {
      if (wspec$weighted) {
        st <- switch(
          nm,
          pp = minimize_pp_blockmodel_dl(
            graph, uniform = FALSE, degree_corrected = degree_corrected,
            degree_prior = degree_prior, B_min = B_min, B_max = B_max,
            refine = refine, weights = wspec$W, types = wspec$types,
            params = params
          ),
          pp_uniform = minimize_pp_blockmodel_dl(
            graph, uniform = TRUE, degree_corrected = degree_corrected,
            degree_prior = degree_prior, B_min = B_min, B_max = B_max,
            refine = refine, weights = wspec$W, types = wspec$types,
            params = params
          ),
          sbm_flat = minimize_weighted_blockmodel_dl(
            graph, weights = wspec$W, types = wspec$types, params = params,
            degree_corrected = degree_corrected, degree_prior = degree_prior,
            B_min = B_min, B_max = B_max, refine = refine
          ),
          sbm_nested = minimize_nested_weighted_blockmodel_dl(
            graph, weights = wspec$W, types = wspec$types, params = params,
            degree_corrected = degree_corrected, degree_prior = degree_prior,
            B_min = B_min, B_max = B_max, refine = refine
          )
        )
      } else {
        st <- switch(
          nm,
          pp = minimize_pp_blockmodel_dl(
            graph, uniform = FALSE, degree_corrected = degree_corrected,
            degree_prior = degree_prior, B_min = B_min, B_max = B_max, refine = refine
          ),
          pp_uniform = minimize_pp_blockmodel_dl(
            graph, uniform = TRUE, degree_corrected = degree_corrected,
            degree_prior = degree_prior, B_min = B_min, B_max = B_max, refine = refine
          ),
          sbm_flat = minimize_blockmodel_dl(
            graph, degree_corrected = degree_corrected,
            degree_prior = degree_prior, B_min = B_min, B_max = B_max, refine = refine
          ),
          sbm_nested = minimize_nested_blockmodel_dl(
            graph, degree_corrected = degree_corrected,
            degree_prior = degree_prior, B_min = B_min, B_max = B_max, refine = refine
          )
        )
      }
      st
    }, cores = cores, graph = graph)
    for (run in seq_len(nruns)) {
      st <- .starts[[run]]
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
    weighted = wspec$weighted,
    types = wspec$types,
    method = method,
    evidence_method = if (method == "posterior") evidence_method else NULL
  ), class = "rgraphtool_assortative_model_comparison")
}

#' @noRd
#' @export
print.rgraphtool_assortative_model_comparison <- function(x, ...) {
  wtxt <- if (isTRUE(x$weighted)) {
    sprintf(", weighted (%s)", paste(x$types, collapse = ", "))
  } else {
    ""
  }
  cat(sprintf("<assortative vs. general SBM comparison (%s, %s%s)>\n",
              if (x$degree_corrected) "DC" else "NDC",
              if (x$method == "posterior") paste0("posterior/", x$evidence_method) else "MAP",
              wtxt))
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
