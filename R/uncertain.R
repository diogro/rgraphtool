# Uncertain network reconstruction, measured networks, extraneous edge
# probabilities, latent Poisson multigraphs, and Bayesian edge prediction.
#
# Copyright (C) 2006-2026 Tiago de Paula Peixoto (graph-tool algorithms)
# Copyright (C) 2026 Diogo Melo (R port)
# LGPL-3+

.MODEL_MEASURED <- 0L
.MODEL_MIXED <- 1L
.MODEL_UNCERTAIN <- 2L
.MODEL_LATENT_MULTI <- 3L

.uncertain_resolve_edge_vec <- function(graph, val, default_val, name) {
  m <- graph_num_edges(graph)
  if (is.null(val)) {
    return(rep(default_val, m))
  }
  if (is.character(val) && length(val) == 1L) {
    v <- graph_edge_property(graph, val)
    return(as.numeric(v))
  }
  if (!is.numeric(val) || length(val) != m || anyNA(val)) {
    stop(sprintf("%s must be NULL, an edge property name, or a numeric vector of length %d",
                 name, m), call. = FALSE)
  }
  as.numeric(val)
}

.uncertain_validate_bs <- function(N, blocks, nested) {
  if (is.null(blocks)) {
    if (nested) return(list(rep(1L, N), 1L))
    return(list(rep(1L, N)))
  }
  if (inherits(blocks, "rgraphtool_nested_block_state")) {
    bs <- lapply(blocks$levels, function(x) as.integer(x))
    if (!nested) return(list(bs[[1L]]))
    return(bs)
  }
  if (inherits(blocks, "rgraphtool_block_state")) {
    b0 <- as.integer(blocks$blocks)
    if (nested) return(list(b0, rep(1L, max(1L, max(b0)))))
    return(list(b0))
  }
  if (is.list(blocks)) {
    bs <- lapply(blocks, function(x) {
      if (!is.numeric(x) || anyNA(x) || any(x < 1)) {
        stop("Each element of blocks must be a valid 1-based integer vector", call. = FALSE)
      }
      as.integer(x)
    })
    if (length(bs) == 0L || length(bs[[1L]]) != N) {
      stop(sprintf("blocks[[1]] must have length N (%d)", N), call. = FALSE)
    }
    if (!nested && length(bs) > 1L) bs <- bs[1L]
    if (nested && length(bs) == 1L) {
      bs[[2L]] <- rep(1L, max(1L, max(bs[[1L]])))
    }
    return(bs)
  }
  if (is.numeric(blocks) && length(blocks) == N && !anyNA(blocks) && all(blocks >= 1)) {
    b0 <- as.integer(blocks)
    if (nested) return(list(b0, rep(1L, max(1L, max(b0)))))
    return(list(b0))
  }
  stop("blocks must be NULL, a length-N integer vector, a hierarchy list, or a block state",
       call. = FALSE)
}

.uncertain_validate_edges_mat <- function(edges, N, name = "edges") {
  if (is.data.frame(edges)) {
    if (ncol(edges) < 2L) stop(sprintf("%s must have at least 2 columns", name), call. = FALSE)
    edges <- as.matrix(edges[, 1:2, drop = FALSE])
  }
  if (!is.matrix(edges) || ncol(edges) != 2L || !is.numeric(edges)) {
    stop(sprintf("%s must be a two-column numeric matrix or data.frame", name), call. = FALSE)
  }
  if (nrow(edges) > 0L && (anyNA(edges) || any(edges < 1) || any(edges > N) ||
                           any(edges != floor(edges)))) {
    stop(sprintf("%s must contain valid 1-based vertex IDs in 1..%d", name, N), call. = FALSE)
  }
  storage.mode(edges) <- "integer"
  edges
}

.uncertain_build_state <- function(model_code, model_name, graph, obs_edges,
                                   obs_v1, obs_v2, def1, def2, hparams,
                                   lp, lq, nested, blocks, latent_edges,
                                   latent_weights, degree_corrected, degree_prior,
                                   self_loops, q_max) {
  graph_pointer(graph)
  N <- graph_num_vertices(graph)
  directed <- isTRUE(graph$directed)
  degree_prior <- match.arg(degree_prior, c("distributed", "uniform"))
  dist_flag <- identical(degree_prior, "distributed")
  nested <- isTRUE(nested)
  degree_corrected <- isTRUE(degree_corrected)
  self_loops <- isTRUE(self_loops)

  bs <- .uncertain_validate_bs(N, blocks, nested)

  if (is.null(latent_edges)) {
    lat_e <- obs_edges
    lat_w <- if (is.null(latent_weights)) rep(1L, nrow(lat_e)) else as.integer(latent_weights)
  } else {
    lat_e <- .uncertain_validate_edges_mat(latent_edges, N, "latent_edges")
    lat_w <- if (is.null(latent_weights)) rep(1L, nrow(lat_e)) else as.integer(latent_weights)
  }
  if (length(lat_w) != nrow(lat_e) || anyNA(lat_w) || any(lat_w < 0L)) {
    stop("latent_weights must be a non-negative integer vector matching latent_edges", call. = FALSE)
  }
  keep <- lat_w > 0L
  lat_e <- lat_e[keep, , drop = FALSE]
  lat_w <- lat_w[keep]

  structure(
    list(
      model_code = as.integer(model_code),
      model = model_name,
      graph = graph,
      n_vertices = as.integer(N),
      directed = directed,
      self_loops = self_loops,
      nested = nested,
      degree_corrected = degree_corrected,
      degree_prior = degree_prior,
      dist = dist_flag,
      q_max = as.numeric(q_max),
      obs_edges = obs_edges,
      obs_v1 = as.numeric(obs_v1),
      obs_v2 = as.numeric(obs_v2),
      def1 = as.numeric(def1),
      def2 = as.numeric(def2),
      hparams = as.numeric(hparams),
      lp = as.numeric(lp),
      lq = as.numeric(lq),
      latent_edges = lat_e,
      latent_weights = as.integer(lat_w),
      bs = bs
    ),
    class = "rgraphtool_uncertain_block_state"
  )
}

#' @noRd
#' @export
measured_block_state <- function(graph, n = NULL, x = NULL,
                                 n_default = 1L, x_default = 0L,
                                 fn_params = list(alpha = 1, beta = 1),
                                 fp_params = list(mu = 1, nu = 1),
                                 lp = NA_real_, lq = NA_real_,
                                 nested = TRUE, blocks = NULL,
                                 latent_edges = NULL, latent_weights = NULL,
                                 degree_corrected = TRUE,
                                 degree_prior = c("distributed", "uniform"),
                                 self_loops = FALSE, q_max = Inf) {
  graph_pointer(graph)
  obs_edges <- graph_edges(graph)
  nv <- .uncertain_resolve_edge_vec(graph, n, 1, "n")
  xv <- .uncertain_resolve_edge_vec(graph, x, 1, "x")
  if (any(nv < 0) || any(xv < 0) || any(xv > nv)) {
    stop("Measured counts must satisfy 0 <= x <= n", call. = FALSE)
  }
  if (n_default < 0 || x_default < 0 || x_default > n_default) {
    stop("Default counts must satisfy 0 <= x_default <= n_default", call. = FALSE)
  }
  alpha <- if (!is.null(fn_params$alpha)) fn_params$alpha else 1
  beta  <- if (!is.null(fn_params$beta))  fn_params$beta  else 1
  mu    <- if (!is.null(fp_params$mu))    fp_params$mu    else 1
  nu    <- if (!is.null(fp_params$nu))    fp_params$nu    else 1
  hparams <- c(alpha, beta, mu, nu)
  if (any(!is.finite(hparams)) || any(hparams <= 0)) {
    stop("Beta hyperparameters (alpha, beta, mu, nu) must be positive finite numbers", call. = FALSE)
  }
  if (!is.na(lp) && (!is.finite(lp) || lp > 0)) {
    stop("lp must be NA or a non-positive log-probability (<= 0)", call. = FALSE)
  }
  if (!is.na(lq) && (!is.finite(lq) || lq > 0)) {
    stop("lq must be NA or a non-positive log-probability (<= 0)", call. = FALSE)
  }
  .uncertain_build_state(
    .MODEL_MEASURED, "measured", graph, obs_edges, nv, xv,
    n_default, x_default, hparams,
    if (is.na(lp)) NaN else lp, if (is.na(lq)) NaN else lq,
    nested, blocks, latent_edges, latent_weights,
    degree_corrected, match.arg(degree_prior), self_loops, q_max
  )
}

#' @noRd
#' @export
mixed_measured_block_state <- function(graph, n = NULL, x = NULL,
                                       n_default = 1L, x_default = 0L,
                                       fn_params = list(alpha = 1, beta = 1),
                                       fp_params = list(mu = 1, nu = 1),
                                       nested = TRUE, blocks = NULL,
                                       latent_edges = NULL, latent_weights = NULL,
                                       degree_corrected = TRUE,
                                       degree_prior = c("distributed", "uniform"),
                                       self_loops = FALSE, q_max = Inf) {
  graph_pointer(graph)
  obs_edges <- graph_edges(graph)
  nv <- .uncertain_resolve_edge_vec(graph, n, 1, "n")
  xv <- .uncertain_resolve_edge_vec(graph, x, 1, "x")
  if (any(nv < 0) || any(xv < 0) || any(xv > nv)) {
    stop("Measured counts must satisfy 0 <= x <= n", call. = FALSE)
  }
  if (n_default < 0 || x_default < 0 || x_default > n_default) {
    stop("Default counts must satisfy 0 <= x_default <= n_default", call. = FALSE)
  }
  alpha <- if (!is.null(fn_params$alpha)) fn_params$alpha else 1
  beta  <- if (!is.null(fn_params$beta))  fn_params$beta  else 1
  mu    <- if (!is.null(fp_params$mu))    fp_params$mu    else 1
  nu    <- if (!is.null(fp_params$nu))    fp_params$nu    else 1
  hparams <- c(alpha, beta, mu, nu)
  if (any(!is.finite(hparams)) || any(hparams <= 0)) {
    stop("Beta hyperparameters (alpha, beta, mu, nu) must be positive finite numbers", call. = FALSE)
  }
  .uncertain_build_state(
    .MODEL_MIXED, "mixed_measured", graph, obs_edges, nv, xv,
    n_default, x_default, hparams, NaN, NaN,
    nested, blocks, latent_edges, latent_weights,
    degree_corrected, match.arg(degree_prior), self_loops, q_max
  )
}

#' @noRd
#' @export
uncertain_block_state <- function(graph, q = NULL, q_default = 0,
                                  nested = TRUE, blocks = NULL,
                                  latent_edges = NULL, latent_weights = NULL,
                                  degree_corrected = TRUE,
                                  degree_prior = c("distributed", "uniform"),
                                  self_loops = FALSE, q_max = Inf) {
  graph_pointer(graph)
  obs_edges <- graph_edges(graph)
  qv <- .uncertain_resolve_edge_vec(graph, q, 0.5, "q")
  if (any(qv < 0) || any(qv > 1)) {
    stop("Edge probabilities q must lie in [0, 1]", call. = FALSE)
  }
  if (!is.numeric(q_default) || length(q_default) != 1L || is.na(q_default) ||
      q_default < 0 || q_default > 1) {
    stop("q_default must be a single numeric value in [0, 1]", call. = FALSE)
  }
  .uncertain_build_state(
    .MODEL_UNCERTAIN, "uncertain", graph, obs_edges, qv, numeric(length(qv)),
    q_default, 0, c(1, 1, 1, 1), NaN, NaN,
    nested, blocks, latent_edges, latent_weights,
    degree_corrected, match.arg(degree_prior), self_loops, q_max
  )
}

#' @noRd
#' @export
latent_multigraph_block_state <- function(graph, nested = TRUE, blocks = NULL,
                                          latent_weights = NULL,
                                          degree_corrected = TRUE,
                                          degree_prior = c("distributed", "uniform"),
                                          self_loops = TRUE, q_max = Inf) {
  graph_pointer(graph)
  obs_edges <- graph_edges(graph)
  m <- nrow(obs_edges)
  .uncertain_build_state(
    .MODEL_LATENT_MULTI, "latent_multigraph", graph, obs_edges,
    numeric(m), numeric(m), 0, 0, c(1, 1, 1, 1), NaN, NaN,
    nested, blocks, obs_edges, latent_weights,
    degree_corrected, match.arg(degree_prior), self_loops, q_max
  )
}

.check_uncertain_state <- function(state) {
  if (!inherits(state, "rgraphtool_uncertain_block_state")) {
    stop("state must be created by measured_block_state(), mixed_measured_block_state(), uncertain_block_state(), or latent_multigraph_block_state()",
         call. = FALSE)
  }
}

#' @noRd
#' @export
uncertain_description_length <- function(state, sbm = TRUE, density = FALSE,
                                         aE = 1, breakdown = FALSE) {
  .check_uncertain_state(state)
  res <- cpp_uncertain_entropy(
    state$model_code, state$n_vertices, state$directed, state$self_loops,
    state$obs_edges, state$obs_v1, state$obs_v2,
    state$def1, state$def2, state$hparams, state$lp, state$lq,
    state$latent_edges, state$latent_weights,
    state$bs, state$degree_corrected, state$dist,
    isTRUE(sbm), isTRUE(density), as.numeric(aE), state$q_max
  )
  if (isTRUE(breakdown)) return(res)
  as.numeric(res$total)
}

#' @noRd
#' @export
measured_error_posterior <- function(state) {
  .check_uncertain_state(state)
  if (!identical(state$model, "measured")) {
    stop("measured_error_posterior() requires a measured_block_state()", call. = FALSE)
  }
  res <- cpp_uncertain_entropy(
    state$model_code, state$n_vertices, state$directed, state$self_loops,
    state$obs_edges, state$obs_v1, state$obs_v2,
    state$def1, state$def2, state$hparams, state$lp, state$lq,
    state$latent_edges, state$latent_weights,
    state$bs, state$degree_corrected, state$dist,
    FALSE, FALSE, 1.0, state$q_max
  )
  pa <- res$p_posterior[1L]; pb <- res$p_posterior[2L]
  qa <- res$q_posterior[1L]; qb <- res$q_posterior[2L]
  list(
    p_alpha = pa, p_beta = pb,
    p_mean = pa / (pa + pb),
    p_sd = sqrt(pa * pb / ((pa + pb)^2 * (pa + pb + 1))),
    q_alpha = qa, q_beta = qb,
    q_mean = qa / (qa + qb),
    q_sd = sqrt(qa * qb / ((qa + qb)^2 * (qa + qb + 1)))
  )
}

#' @noRd
#' @export
edge_probabilities <- function(state, edges, log = TRUE, epsilon = 1e-8) {
  .check_uncertain_state(state)
  qmat <- .uncertain_validate_edges_mat(edges, state$n_vertices, "edges")
  lp <- cpp_uncertain_edge_probs(
    state$model_code, state$n_vertices, state$directed, state$self_loops,
    state$obs_edges, state$obs_v1, state$obs_v2,
    state$def1, state$def2, state$hparams, state$lp, state$lq,
    state$latent_edges, state$latent_weights,
    state$bs, state$degree_corrected, state$dist,
    qmat, as.numeric(epsilon), state$q_max
  )
  if (isTRUE(log)) lp else exp(lp)
}

.uncertain_init_sbm_bs <- function(state) {
  N <- state$n_vertices
  if (N <= 1L || nrow(state$latent_edges) == 0L) return(state$bs)
  # Only auto-initialize if currently a trivial single-block partition
  if (length(unique(state$bs[[1L]])) > 1L) return(state$bs)
  g_lat <- graph_create(N, directed = state$directed)
  if (nrow(state$latent_edges) > 0L) {
    # Expand multiplicity (cap at 10 for fast initialization)
    w <- pmin(state$latent_weights, 10L)
    idx <- rep.int(seq_len(nrow(state$latent_edges)), w)
    graph_add_edges(g_lat, state$latent_edges[idx, , drop = FALSE])
  }
  if (state$nested) {
    fit <- minimize_nested_blockmodel_dl(
      g_lat, degree_corrected = state$degree_corrected,
      degree_prior = state$degree_prior
    )
    lapply(fit$levels, as.integer)
  } else {
    fit <- minimize_blockmodel_dl(
      g_lat, degree_corrected = state$degree_corrected,
      degree_prior = state$degree_prior
    )
    list(as.integer(fit$blocks))
  }
}

#' @noRd
#' @export
reconstruct_network <- function(state, nsamples = 200L, burnin = 100L, thin = 1L,
                                beta = 1, pedges = 0.5, ph = 0.1, hstep = 1,
                                c = 1, d = 0.01, init_sbm = TRUE,
                                collect_partitions = FALSE, threshold = 0.5) {
  .check_uncertain_state(state)
  nsamples <- as.integer(nsamples)
  burnin <- as.integer(burnin)
  thin <- as.integer(thin)
  if (is.na(nsamples) || nsamples < 1L) stop("nsamples must be >= 1", call. = FALSE)
  if (is.na(burnin) || burnin < 0L) stop("burnin must be >= 0", call. = FALSE)
  if (is.na(thin) || thin < 1L) stop("thin must be >= 1", call. = FALSE)

  bs_start <- if (isTRUE(init_sbm)) .uncertain_init_sbm_bs(state) else state$bs
  res <- cpp_uncertain_reconstruct(
    state$model_code, state$n_vertices, state$directed, state$self_loops,
    state$obs_edges, state$obs_v1, state$obs_v2,
    state$def1, state$def2, state$hparams, state$lp, state$lq,
    state$latent_edges, state$latent_weights,
    bs_start, state$degree_corrected, state$dist,
    nsamples, burnin, thin, as.numeric(beta),
    as.numeric(pedges), as.numeric(ph), as.numeric(hstep),
    as.numeric(c), as.numeric(d), isTRUE(collect_partitions), state$q_max
  )

  edge_df <- data.frame(
    from = res$from,
    to = res$to,
    prob = res$prob,
    mean_weight = res$mean_weight,
    stringsAsFactors = FALSE
  )
  if (nrow(edge_df) > 0L) {
    ord <- order(-edge_df$prob, edge_df$from, edge_df$to)
    edge_df <- edge_df[ord, , drop = FALSE]
    rownames(edge_df) <- NULL
  }

  g_cons <- graph_create(state$n_vertices, directed = state$directed)
  keep <- edge_df$prob >= threshold
  if (any(keep)) {
    sub_df <- edge_df[keep, , drop = FALSE]
    emat <- cbind(as.integer(sub_df$from), as.integer(sub_df$to))
    graph_add_edges(g_cons, emat, properties = list(
      prob = sub_df$prob,
      mean_weight = sub_df$mean_weight
    ))
  }

  new_state <- state
  new_state$latent_edges <- res$final_edges
  new_state$latent_weights <- as.integer(res$final_weights)
  new_state$bs <- res$final_bs
  new_state$hparams <- res$final_hparams

  out <- list(
    model = state$model,
    edges = edge_df,
    consensus_graph = g_cons,
    threshold = threshold,
    state = new_state,
    p_posterior = res$p_posterior,
    q_posterior = res$q_posterior,
    p_samples = res$p_samples,
    q_samples = res$q_samples,
    final_hparams = res$final_hparams,
    entropies = res$entropies,
    entropy = res$entropy,
    dS = res$dS,
    nattempts = res$nattempts,
    nmoves = res$nmoves
  )
  if (isTRUE(collect_partitions)) out$partitions <- res$partitions
  structure(out, class = "rgraphtool_reconstruction")
}

#' @noRd
#' @export
latent_multigraph <- function(graph, epsilon = 1e-8, max_niter = 0L) {
  ptr <- graph_pointer(graph)
  if (!is.numeric(epsilon) || length(epsilon) != 1L || is.na(epsilon) || epsilon <= 0) {
    stop("epsilon must be a positive number", call. = FALSE)
  }
  max_niter <- as.integer(max_niter)
  if (is.na(max_niter) || max_niter < 0L) stop("max_niter must be >= 0", call. = FALSE)
  res <- cpp_latent_multigraph_em(ptr, as.numeric(epsilon), max_niter)
  edf <- data.frame(
    from = res$edges[, 1L],
    to = res$edges[, 2L],
    weight = res$weights,
    is_self_loop = res$edges[, 1L] == res$edges[, 2L],
    stringsAsFactors = FALSE
  )
  structure(
    list(
      edges = edf,
      theta_out = res$theta_out,
      theta_in = res$theta_in,
      niter = res$niter,
      delta = res$delta,
      directed = isTRUE(graph$directed)
    ),
    class = "rgraphtool_latent_multigraph"
  )
}

#' @noRd
#' @export
predict_edges <- function(x, candidates, nsamples = 100L, burnin = 50L, thin = 1L,
                          nested = TRUE, degree_corrected = TRUE,
                          degree_prior = c("distributed", "uniform"),
                          self_loops = FALSE, epsilon = 1e-8) {
  degree_prior <- match.arg(degree_prior)
  if (inherits(x, "rgraphtool_uncertain_block_state")) {
    N <- x$n_vertices
    cand_mat <- .uncertain_validate_edges_mat(candidates, N, "candidates")
    rec <- reconstruct_network(
      x, nsamples = nsamples, burnin = burnin, thin = thin,
      init_sbm = TRUE, collect_partitions = FALSE
    )
    cond_lp <- edge_probabilities(rec$state, cand_mat, log = TRUE, epsilon = epsilon)
    key_fn <- function(u, v, dir) {
      if (!dir) {
        a <- pmin(u, v); b <- pmax(u, v)
      } else {
        a <- u; b <- v
      }
      paste0(a, ":", b)
    }
    m_keys <- key_fn(rec$edges$from, rec$edges$to, x$directed)
    c_keys <- key_fn(cand_mat[, 1L], cand_mat[, 2L], x$directed)
    idx <- match(c_keys, m_keys)
    marg_prob <- ifelse(is.na(idx), 0, rec$edges$prob[idx])
    res_df <- data.frame(
      from = cand_mat[, 1L],
      to = cand_mat[, 2L],
      prob = marg_prob,
      conditional_prob = exp(cond_lp),
      log_prob = cond_lp,
      stringsAsFactors = FALSE
    )
    return(structure(
      list(predictions = res_df, model = x$model, nsamples = as.integer(nsamples),
           reconstruction = rec),
      class = "rgraphtool_edge_prediction"
    ))
  }

  # Graph or fitted SBM state: pure SBM posterior edge prediction
  if (inherits(x, "rgraphtool_nested_block_state")) {
    graph <- x$graph
    nested <- TRUE
    degree_corrected <- isTRUE(x$degree_corrected)
    degree_prior <- x$degree_prior
    fit <- x
  } else if (inherits(x, "rgraphtool_block_state")) {
    graph <- x$graph
    nested <- FALSE
    degree_corrected <- isTRUE(x$degree_corrected)
    degree_prior <- x$degree_prior
    fit <- x
  } else {
    graph_pointer(x)
    graph <- x
    fit <- if (isTRUE(nested)) {
      minimize_nested_blockmodel_dl(graph, degree_corrected = degree_corrected,
                                    degree_prior = degree_prior)
    } else {
      minimize_blockmodel_dl(graph, degree_corrected = degree_corrected,
                             degree_prior = degree_prior)
    }
  }

  N <- graph_num_vertices(graph)
  directed <- isTRUE(graph$directed)
  cand_mat <- .uncertain_validate_edges_mat(candidates, N, "candidates")
  K <- nrow(cand_mat)
  obs_e <- graph_edges(graph)

  nsamples <- as.integer(nsamples)
  burnin <- as.integer(burnin)
  thin <- as.integer(thin)
  if (is.na(nsamples) || nsamples < 1L) stop("nsamples must be >= 1", call. = FALSE)

  post <- sample_posterior(fit, nsamples = nsamples, burnin = burnin, thin = thin)
  parts <- post$samples

  # Under MODEL_UNCERTAIN with q = 0.5 everywhere, measurement log-odds are identically 0,
  # so cpp_uncertain_edge_probs computes exact pure-SBM log P(w_uv >= 1 | G_{\uv}, b).
  empty_obs <- matrix(integer(0), nrow = 0L, ncol = 2L)
  obs_w <- rep(1L, nrow(obs_e))
  dist_flag <- identical(degree_prior, "distributed")

  prob_sum <- numeric(K)
  lprob_mat <- matrix(NA_real_, nrow = K, ncol = nsamples)
  for (s in seq_len(nsamples)) {
    bs_s <- if (isTRUE(nested)) {
      lapply(parts[[s]], as.integer)
    } else {
      list(as.integer(parts[[s]]))
    }
    lp_s <- cpp_uncertain_edge_probs(
      .MODEL_UNCERTAIN, N, directed, isTRUE(self_loops),
      empty_obs, numeric(0), numeric(0),
      0.5, 0.0, c(1, 1, 1, 1), NaN, NaN,
      obs_e, obs_w, bs_s,
      isTRUE(degree_corrected), dist_flag,
      cand_mat, as.numeric(epsilon), Inf
    )
    lprob_mat[, s] <- lp_s
    prob_sum <- prob_sum + exp(lp_s)
  }
  mean_prob <- prob_sum / nsamples
  # Exact log-mean-exp for numerical stability
  max_lp <- apply(lprob_mat, 1L, max)
  log_mean_prob <- ifelse(
    is.finite(max_lp),
    max_lp + log(rowMeans(exp(lprob_mat - max_lp))),
    -Inf
  )

  # Mark which candidates are currently observed in graph
  ek <- function(u, v) {
    if (!directed) paste0(pmin(u, v), ":", pmax(u, v)) else paste0(u, ":", v)
  }
  obs_set <- unique(ek(obs_e[, 1L], obs_e[, 2L]))
  is_obs <- ek(cand_mat[, 1L], cand_mat[, 2L]) %in% obs_set

  res_df <- data.frame(
    from = cand_mat[, 1L],
    to = cand_mat[, 2L],
    observed = is_obs,
    prob = mean_prob,
    log_prob = log_mean_prob,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      predictions = res_df,
      model = if (isTRUE(nested)) "nested_sbm" else "flat_sbm",
      nsamples = nsamples,
      fit = fit
    ),
    class = "rgraphtool_edge_prediction"
  )
}

#' @noRd
#' @export
evaluate_edge_prediction <- function(scores, labels, threshold = 0.5) {
  if (inherits(scores, "rgraphtool_edge_prediction")) {
    scores <- scores$predictions$prob
  }
  if (!is.numeric(scores) || anyNA(scores)) {
    stop("scores must be a numeric vector without NAs", call. = FALSE)
  }
  if (is.logical(labels)) labels <- as.integer(labels)
  if (!is.numeric(labels) || length(labels) != length(scores) ||
      anyNA(labels) || !all(labels %in% c(0, 1))) {
    stop("labels must be a binary (0/1 or logical) vector matching scores", call. = FALSE)
  }
  n <- length(scores)
  n_pos <- sum(labels == 1)
  n_neg <- sum(labels == 0)
  if (n_pos == 0 || n_neg == 0) {
    stop("labels must contain at least one positive (1) and one negative (0)", call. = FALSE)
  }

  # Exact ROC-AUC via Mann-Whitney U with exact tie averaging
  ranks <- rank(scores, ties.method = "average")
  roc_auc <- (sum(ranks[labels == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)

  # Build exact ROC & PR curves across unique score thresholds
  ord <- order(-scores)
  s_ord <- scores[ord]
  y_ord <- labels[ord]
  tp_cum <- cumsum(y_ord == 1)
  fp_cum <- cumsum(y_ord == 0)
  # Keep last index of each tied score group
  last_in_tie <- c(which(diff(s_ord) != 0), n)
  tp_u <- tp_cum[last_in_tie]
  fp_u <- fp_cum[last_in_tie]
  thr_u <- s_ord[last_in_tie]

  tpr <- c(0, tp_u / n_pos)
  fpr <- c(0, fp_u / n_neg)
  roc_df <- data.frame(threshold = c(Inf, thr_u), fpr = fpr, tpr = tpr)

  prec_u <- tp_u / (tp_u + fp_u)
  rec_u <- tp_u / n_pos
  pr_df <- data.frame(threshold = c(Inf, thr_u), recall = c(0, rec_u),
                      precision = c(1, prec_u))
  # Step-function / trapezoidal PR-AUC (average precision)
  drec <- diff(c(0, rec_u))
  pr_auc <- sum(drec * prec_u)

  # Thresholded confusion matrix metrics
  pred <- as.integer(scores >= threshold)
  tp <- sum(pred == 1L & labels == 1)
  fp <- sum(pred == 1L & labels == 0)
  tn <- sum(pred == 0L & labels == 0)
  fn <- sum(pred == 0L & labels == 1)
  precision <- if (tp + fp > 0) tp / (tp + fp) else NA_real_
  recall <- tp / (tp + fn)
  f1 <- if (!is.na(precision) && (precision + recall) > 0) {
    2 * precision * recall / (precision + recall)
  } else NA_real_
  accuracy <- (tp + tn) / n
  brier <- mean((scores - labels)^2)

  structure(
    list(
      roc_auc = roc_auc,
      pr_auc = pr_auc,
      threshold = threshold,
      precision = precision,
      recall = recall,
      f1 = f1,
      accuracy = accuracy,
      brier = brier,
      confusion = c(tp = tp, fp = fp, tn = tn, fn = fn),
      roc = roc_df,
      pr = pr_df
    ),
    class = "rgraphtool_edge_evaluation"
  )
}

#' @noRd
#' @export
print.rgraphtool_uncertain_block_state <- function(x, ...) {
  cat(sprintf("<rgraphtool_uncertain_block_state (%s, %s, %s)>\n",
              x$model,
              if (x$directed) "directed" else "undirected",
              if (x$nested) "nested" else "flat"))
  cat(sprintf("  Vertices: %d | Observed dyads: %d | Latent edges: %d (total weight %d)\n",
              x$n_vertices, nrow(x$obs_edges), nrow(x$latent_edges),
              sum(x$latent_weights)))
  cat(sprintf("  Blocks (level 1): %d | Degree-corrected: %s\n",
              length(unique(x$bs[[1L]])), if (x$degree_corrected) "TRUE" else "FALSE"))
  invisible(x)
}

#' @noRd
#' @export
print.rgraphtool_reconstruction <- function(x, ...) {
  cat(sprintf("<rgraphtool_reconstruction (%s)>\n", x$model))
  cat(sprintf("  Marginal edges observed: %d | Consensus edges (p >= %.2f): %d\n",
              nrow(x$edges), x$threshold, graph_num_edges(x$consensus_graph)))
  cat(sprintf("  Final description length: %.4f | Moves: %.0f / %.0f\n",
              x$entropy, x$nmoves, x$nattempts))
  if (identical(x$model, "measured") && length(x$p_samples) > 0L) {
    cat(sprintf("  Posterior mean FN rate (p): %.4f | FP rate (q): %.4f\n",
                mean(x$p_samples), mean(x$q_samples)))
  }
  invisible(x)
}

#' @noRd
#' @export
print.rgraphtool_latent_multigraph <- function(x, ...) {
  n_non <- sum(!x$edges$is_self_loop)
  n_sl <- sum(x$edges$is_self_loop & x$edges$weight > 0)
  cat(sprintf("<rgraphtool_latent_multigraph (%s)>\n",
              if (x$directed) "directed" else "undirected"))
  cat(sprintf("  Simple edges: %d | Nonzero self-loops: %d | Total expected weight: %.3f\n",
              n_non, n_sl, sum(x$edges$weight)))
  cat(sprintf("  EM iterations: %d | Final delta: %.3e\n", x$niter, x$delta))
  invisible(x)
}

#' @noRd
#' @export
print.rgraphtool_edge_prediction <- function(x, ...) {
  cat(sprintf("<rgraphtool_edge_prediction (%s, %d posterior samples)>\n",
              x$model, x$nsamples))
  cat(sprintf("  Candidate dyads: %d | Mean predicted prob: %.4f\n",
              nrow(x$predictions), mean(x$predictions$prob)))
  invisible(x)
}

#' @noRd
#' @export
print.rgraphtool_edge_evaluation <- function(x, ...) {
  cat("<rgraphtool_edge_evaluation>\n")
  cat(sprintf("  ROC-AUC: %.4f | PR-AUC: %.4f | Brier: %.4f\n",
              x$roc_auc, x$pr_auc, x$brier))
  cat(sprintf("  At threshold %.2f -> Precision: %.4f | Recall: %.4f | F1: %.4f | Acc: %.4f\n",
              x$threshold, x$precision, x$recall, x$f1, x$accuracy))
  invisible(x)
}
