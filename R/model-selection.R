logsumexp_vec <- function(x) {
  x <- x[is.finite(x) | x == Inf]
  if (!length(x)) return(-Inf)
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

state_description_length <- function(state) {
  if (inherits(state, "graphtool_block_state")) {
    return(block_description_length(state))
  }
  if (inherits(state, "graphtool_nested_block_state")) {
    return(nested_description_length(state))
  }
  if (inherits(state, "graphtool_overlap_block_state")) {
    return(overlap_description_length(state))
  }
  if (inherits(state, "graphtool_nested_overlap_block_state")) {
    return(nested_overlap_description_length(state))
  }
  if (inherits(state, "graphtool_weighted_block_state")) {
    return(weighted_description_length(state))
  }
  if (inherits(state, "graphtool_nested_weighted_block_state")) {
    return(nested_weighted_description_length(state))
  }
  stop("state must be a block, nested, overlap, nested overlap, weighted, or nested weighted block state",
       call. = FALSE)
}

state_num_groups <- function(state) {
  if (inherits(state, "graphtool_block_state")) {
    return(length(unique(block_membership(state))))
  }
  if (inherits(state, "graphtool_nested_block_state")) {
    return(length(unique(nested_levels(state)[[1L]])))
  }
  if (inherits(state, "graphtool_overlap_block_state") ||
      inherits(state, "graphtool_nested_overlap_block_state")) {
    return(as.integer(overlap_groups(state)$B))
  }
  if (inherits(state, "graphtool_weighted_block_state")) {
    return(as.integer(state$B))
  }
  if (inherits(state, "graphtool_nested_weighted_block_state")) {
    return(as.integer(state$B[1L]))
  }
  NA_integer_
}

state_num_levels <- function(state) {
  if (inherits(state, "graphtool_block_state") ||
      inherits(state, "graphtool_overlap_block_state") ||
      inherits(state, "graphtool_weighted_block_state")) {
    return(1L)
  }
  if (inherits(state, "graphtool_nested_block_state")) {
    return(length(nested_levels(state)))
  }
  if (inherits(state, "graphtool_nested_overlap_block_state")) {
    return(length(state$upper_levels) + 1L)
  }
  if (inherits(state, "graphtool_nested_weighted_block_state")) {
    return(length(state$levels))
  }
  NA_integer_
}

model_evidence <- function(state, nsamples = 200L, thin = 5L, burnin = NULL,
                           beta = 1, mle = TRUE,
                           method = c("variational", "importance"),
                           samples = NULL, entropies = NULL,
                           multiflip = TRUE,
                           c = 0.5, d = 0.01, B_min = 1L, B_max = Inf) {
  method <- match.arg(method)
  if (!is.logical(multiflip) || length(multiflip) != 1L || is.na(multiflip)) {
    stop("multiflip must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(mle) || length(mle) != 1L || is.na(mle)) {
    stop("mle must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.numeric(nsamples) || length(nsamples) != 1L || is.na(nsamples) ||
      nsamples < 2 || nsamples != floor(nsamples)) {
    stop("nsamples must be an integer >= 2", call. = FALSE)
  }
  if (!is.numeric(thin) || length(thin) != 1L || is.na(thin) ||
      thin < 1 || thin != floor(thin)) {
    stop("thin must be a positive integer", call. = FALSE)
  }
  if (!is.numeric(beta) || length(beta) != 1L || is.na(beta) || beta <= 0) {
    stop("beta must be a positive number", call. = FALSE)
  }

  is_flat <- inherits(state, "graphtool_block_state") ||
    inherits(state, "graphtool_weighted_block_state")
  is_nested <- inherits(state, "graphtool_nested_block_state") ||
    inherits(state, "graphtool_nested_weighted_block_state")
  is_ov_flat <- inherits(state, "graphtool_overlap_block_state")
  is_ov_nested <- inherits(state, "graphtool_nested_overlap_block_state")
  if (!is_flat && !is_nested && !is_ov_flat && !is_ov_nested) {
    stop("state must be a block, nested, overlap, nested overlap, weighted, or nested weighted block state",
         call. = FALSE)
  }

  init_dl <- state_description_length(state)
  final_state <- state
  is_nested_samples <- is_nested

  if (is.null(samples) || is.null(entropies)) {
    if (is_flat || is_nested) {
      mcmc_args <- list(beta = beta, c = c, d = d, B_min = B_min, B_max = B_max)
      eq_args <- list(wait = 20L, nbreaks = 2L, max_niter = 200L)
      post <- sample_posterior(state, nsamples = as.integer(nsamples),
                               thin = as.integer(thin), burnin = burnin,
                               multiflip = multiflip, mcmc_args = mcmc_args,
                               equilibrate_args = eq_args)
      samples <- post$samples
      entropies <- post$entropies
      final_state <- post$state
      is_nested_samples <- isTRUE(post$nested)
    } else {
      # Flat or nested overlap state
      g <- state$graph
      if (graph_num_edges(g) == 0L) {
        return(structure(list(
          log_evidence = 0,
          neg_log_evidence = 0,
          variational_log_evidence = 0,
          importance_log_evidence = 0,
          mean_dl = 0,
          sd_dl = 0,
          map_dl = 0,
          posterior_entropy = 0,
          method = method,
          nsamples = 0L,
          cluster = NULL,
          samples = list(),
          entropies = numeric(0),
          state = state
        ), class = "graphtool_model_evidence"))
      }
      burn_n <- if (is.null(burnin)) 30L else as.integer(floor(burnin))
      cur_state <- if (burn_n > 0L) {
        overlap_mcmc_sweep(state, beta = beta, niter = burn_n, c = c, d = d,
                           B_min = B_min, B_max = B_max)$state
      } else {
        state
      }
      post <- overlap_mcmc_sweep(cur_state, beta = beta,
                                 niter = as.integer(nsamples) * as.integer(thin),
                                 c = c, d = d, thin = as.integer(thin),
                                 collect = TRUE, B_min = B_min, B_max = B_max)
      samples <- post$majority_samples
      entropies <- post$entropies
      final_state <- post$state
      is_nested_samples <- FALSE
    }
  } else {
    if (length(samples) != length(entropies) || length(samples) < 2L) {
      stop("samples and entropies must have the same length (>= 2)", call. = FALSE)
    }
    is_nested_samples <- is.list(samples[[1L]])
  }

  entropies <- as.numeric(entropies)
  M <- length(samples)
  mc <- mode_cluster_state(samples, nested = is_nested_samples)
  H_post <- mode_entropy(mc, mle = mle)
  mean_dl <- mean(entropies)
  sd_dl <- stats::sd(entropies)
  map_dl <- min(init_dl, min(entropies))

  # Variational / graph-tool cookbook estimate: log P(A) = -<S> + H(b)
  var_log_ev <- -mean_dl + H_post

  # Peixoto (2021) Eq. 29 mode-mixture importance estimate:
  # -log P(A) = log((1/M) * sum_m exp(S_m + log q(b_m)))
  # Use each sample's exact aligned copy in its assigned mode to avoid any label-tie -Inf
  mode_pos <- integer(mc$K)
  log_q <- numeric(M)
  for (m in seq_len(M)) {
    k_own <- mc$clusters[m]
    mode_pos[k_own] <- mode_pos[k_own] + 1L
    aligned_own <- mc$modes[[k_own]]$partitions[[mode_pos[k_own]]]
    lps <- numeric(mc$K)
    for (k in seq_len(mc$K)) {
      wk <- mc$weights[k]
      if (wk <= 0) {
        lps[k] <- -Inf
      } else if (k == k_own) {
        lps[k] <- log(wk) + mode_lprob(mc$modes[[k]], aligned_own,
                                       mle = mle, align = FALSE)
      } else {
        lps[k] <- log(wk) + mode_lprob(mc$modes[[k]], samples[[m]],
                                       mle = mle, align = TRUE)
      }
    }
    log_q[m] <- logsumexp_vec(lps)
  }
  imp_neg_log_ev <- logsumexp_vec(entropies + log_q) - log(M)
  imp_log_ev <- -imp_neg_log_ev

  log_ev <- if (method == "variational") var_log_ev else imp_log_ev

  structure(list(
    log_evidence = log_ev,
    neg_log_evidence = -log_ev,
    variational_log_evidence = var_log_ev,
    importance_log_evidence = imp_log_ev,
    mean_dl = mean_dl,
    sd_dl = sd_dl,
    map_dl = map_dl,
    posterior_entropy = H_post,
    method = method,
    nsamples = M,
    cluster = mc,
    samples = samples,
    entropies = entropies,
    state = final_state
  ), class = "graphtool_model_evidence")
}

print.graphtool_model_evidence <- function(x, ...) {
  cat(sprintf("<model evidence (%s, M = %d samples)>\n", x$method, x$nsamples))
  cat(sprintf("  log P(A)          : %.4f nats  (-log P(A) = %.4f)\n",
              x$log_evidence, x$neg_log_evidence))
  cat(sprintf("  mean DL <Sigma>   : %.4f +/- %.4f nats  (MAP = %.4f)\n",
              x$mean_dl, x$sd_dl, x$map_dl))
  cat(sprintf("  posterior entropy : %.4f nats  (K = %d modes, Be = %.2f)\n",
              x$posterior_entropy,
              if (is.null(x$cluster)) 0L else x$cluster$K,
              if (is.null(x$cluster)) 0 else x$cluster$Be))
  cat(sprintf("  [variational = %.4f, importance = %.4f]\n",
              x$variational_log_evidence, x$importance_log_evidence))
  invisible(x)
}

model_class_specs <- function() {
  list(
    flat_ndc           = list(nested = FALSE, overlap = FALSE, degree_corrected = FALSE),
    flat_dc            = list(nested = FALSE, overlap = FALSE, degree_corrected = TRUE),
    nested_ndc         = list(nested = TRUE,  overlap = FALSE, degree_corrected = FALSE),
    nested_dc          = list(nested = TRUE,  overlap = FALSE, degree_corrected = TRUE),
    overlap_flat_ndc   = list(nested = FALSE, overlap = TRUE,  degree_corrected = FALSE),
    overlap_flat_dc    = list(nested = FALSE, overlap = TRUE,  degree_corrected = TRUE),
    overlap_nested_ndc = list(nested = TRUE,  overlap = TRUE,  degree_corrected = FALSE),
    overlap_nested_dc  = list(nested = TRUE,  overlap = TRUE,  degree_corrected = TRUE)
  )
}

fit_single_model_class <- function(graph, spec, degree_prior, nruns,
                                   B_min, B_max, refine) {
  best_state <- NULL
  best_dl <- Inf
  for (run in seq_len(nruns)) {
    st <- if (!spec$nested && !spec$overlap) {
      minimize_blockmodel_dl(graph, degree_corrected = spec$degree_corrected,
                             degree_prior = degree_prior, B_min = B_min,
                             B_max = B_max, refine = refine)
    } else if (spec$nested && !spec$overlap) {
      minimize_nested_blockmodel_dl(graph, degree_corrected = spec$degree_corrected,
                                    degree_prior = degree_prior, B_min = B_min,
                                    B_max = B_max, refine = refine)
    } else if (!spec$nested && spec$overlap) {
      minimize_overlap_blockmodel_dl(graph, degree_corrected = spec$degree_corrected,
                                     degree_prior = degree_prior, B_min = B_min,
                                     B_max = B_max, refine = refine)
    } else {
      minimize_nested_overlap_blockmodel_dl(graph, degree_corrected = spec$degree_corrected,
                                            degree_prior = degree_prior,
                                            B_min = B_min, B_max = B_max,
                                            refine = refine)
    }
    dl <- state_description_length(st)
    if (dl < best_dl) {
      best_dl <- dl
      best_state <- st
    }
  }
  best_state
}

compare_model_classes <- function(graph,
                                  models = c("flat_ndc", "flat_dc",
                                             "nested_ndc", "nested_dc",
                                             "overlap_flat_ndc", "overlap_flat_dc",
                                             "overlap_nested_ndc", "overlap_nested_dc"),
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
  if (!is.character(models) || !length(models) || anyNA(models)) {
    stop("models must be a non-empty character vector of model class names",
         call. = FALSE)
  }
  models <- unique(models)
  specs <- model_class_specs()
  bad <- setdiff(models, names(specs))
  if (length(bad)) {
    stop("unknown model class(es): ", paste(bad, collapse = ", "),
         ". Valid choices are: ", paste(names(specs), collapse = ", "),
         call. = FALSE)
  }
  if (!is.numeric(nruns) || length(nruns) != 1L || is.na(nruns) ||
      nruns < 1 || nruns != floor(nruns)) {
    stop("nruns must be a positive integer", call. = FALSE)
  }
  nruns <- as.integer(nruns)

  states <- stats::setNames(vector("list", length(models)), models)
  evidences <- if (method == "posterior") {
    stats::setNames(vector("list", length(models)), models)
  } else {
    NULL
  }

  rows <- vector("list", length(models))
  for (i in seq_along(models)) {
    nm <- models[i]
    sp <- specs[[nm]]
    st <- fit_single_model_class(graph, sp, degree_prior, nruns,
                                 B_min, B_max, refine)
    states[[nm]] <- st
    dl <- state_description_length(st)
    if (method == "posterior") {
      ev <- model_evidence(st, nsamples = nsamples, thin = thin,
                           burnin = burnin, mle = mle,
                           method = evidence_method,
                           B_min = B_min, B_max = B_max)
      evidences[[nm]] <- ev
      neg_log_ev <- ev$neg_log_evidence
      post_h <- ev$posterior_entropy
    } else {
      neg_log_ev <- dl
      post_h <- NA_real_
    }
    rows[[i]] <- data.frame(
      model = nm,
      nested = sp$nested,
      overlap = sp$overlap,
      degree_corrected = sp$degree_corrected,
      B = state_num_groups(st),
      levels = state_num_levels(st),
      description_length = dl,
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
  lbf <- matrix(0, nrow = k, ncol = k,
                dimnames = list(tab$model, tab$model))
  for (i in seq_len(k)) {
    for (j in seq_len(k)) {
      # log P(A | M_i) - log P(A | M_j) = (-log P(A | M_j)) - (-log P(A | M_i))
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
    method = method,
    evidence_method = if (method == "posterior") evidence_method else NULL,
    degree_prior = degree_prior
  ), class = "graphtool_model_class_comparison")
}

print.graphtool_model_class_comparison <- function(x, ...) {
  cat(sprintf("<SBM model class comparison (%s, %s degree prior)>\n",
              if (x$method == "posterior") paste0("posterior/", x$evidence_method) else "MAP description length",
              x$degree_prior))
  df <- x$summary
  show_cols <- if (x$method == "posterior") {
    c("model", "B", "levels", "description_length", "posterior_entropy",
      "neg_log_evidence", "delta", "posterior_weight")
  } else {
    c("model", "B", "levels", "description_length", "delta", "posterior_weight")
  }
  out <- df[, show_cols, drop = FALSE]
  for (col in intersect(c("description_length", "posterior_entropy",
                          "neg_log_evidence", "delta"), names(out))) {
    out[[col]] <- round(out[[col]], 3)
  }
  out$posterior_weight <- signif(out$posterior_weight, 3)
  print(out, row.names = FALSE)
  cat(sprintf("Best model: %s\n", x$best_model))
  invisible(x)
}
