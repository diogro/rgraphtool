is_nested_state <- function(state) {
  inherits(state, "rgraphtool_nested_block_state") ||
    inherits(state, "rgraphtool_nested_weighted_block_state") ||
    inherits(state, "rgraphtool_nested_layered_block_state") ||
    inherits(state, "rgraphtool_nested_ranked_block_state")
}

is_weighted_state <- function(state) {
  inherits(state, "rgraphtool_weighted_block_state") ||
    inherits(state, "rgraphtool_nested_weighted_block_state")
}

is_layered_state <- function(state) {
  inherits(state, "rgraphtool_layered_block_state") ||
    inherits(state, "rgraphtool_nested_layered_block_state")
}

is_pp_state <- function(state) {
  inherits(state, "rgraphtool_pp_block_state")
}

is_ranked_state <- function(state) {
  inherits(state, "rgraphtool_ranked_block_state") ||
    inherits(state, "rgraphtool_nested_ranked_block_state")
}

is_bipartite_state <- function(state) {
  inherits(state, "rgraphtool_bipartite_block_state")
}

validate_mcmc_state <- function(state) {
  if (inherits(state, "rgraphtool_bipartite_block_state")) {
    validate_bipartite_state(state)
  } else if (inherits(state, "rgraphtool_nested_weighted_block_state")) {
    validate_nested_weighted_state(state)
  } else if (inherits(state, "rgraphtool_weighted_block_state")) {
    validate_weighted_state(state)
  } else if (inherits(state, "rgraphtool_nested_layered_block_state")) {
    validate_nested_layered_state(state)
  } else if (inherits(state, "rgraphtool_layered_block_state")) {
    validate_layered_state(state)
  } else if (inherits(state, "rgraphtool_pp_block_state")) {
    validate_pp_state(state)
  } else if (inherits(state, "rgraphtool_nested_ranked_block_state")) {
    validate_nested_ranked_state(state)
  } else if (inherits(state, "rgraphtool_ranked_block_state")) {
    validate_ranked_state(state)
  } else if (inherits(state, "rgraphtool_nested_block_state")) {
    validate_nested_state(state)
  } else {
    validate_block_state(state)
    if (is.null(state$graph) ||
        graph_num_edges(state$graph) != state$statistics$num_edges ||
        graph_num_vertices(state$graph) != length(state$blocks)) {
      stop("the graph has changed since the block state was created", call. = FALSE)
    }
  }
}

build_mcmc_opts <- function(nested, beta = 1, niter = 1L, c = NULL, d = 0.01,
                            multiflip = TRUE, B_min = 1L, B_max = Inf,
                            psingle = NULL, psplit = 1, pmerge = 1,
                            pmergesplit = 1, gibbs_sweeps = 10L,
                            ls_shuffle = TRUE) {
  if (!is.numeric(beta) || length(beta) != 1L || is.na(beta) || beta < 0) {
    stop("beta must be a non-negative number (or Inf)", call. = FALSE)
  }
  if (!is.numeric(niter) || length(niter) != 1L || is.na(niter) || niter < 0 ||
      niter != floor(niter)) {
    stop("niter must be a non-negative integer", call. = FALSE)
  }
  if (!is.numeric(d) || length(d) != 1L || is.na(d) || d < 0 || d >= 1) {
    stop("d must be in [0, 1)", call. = FALSE)
  }
  if (!is.logical(multiflip) || length(multiflip) != 1L || is.na(multiflip)) {
    stop("multiflip must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.numeric(B_min) || length(B_min) != 1L || is.na(B_min) || B_min < 1 ||
      !is.numeric(B_max) || length(B_max) != 1L || is.na(B_max) || B_max < B_min) {
    stop("B_min and B_max must satisfy 1 <= B_min <= B_max", call. = FALSE)
  }
  opts <- list(
    beta = as.numeric(beta),
    d = as.numeric(d),
    niter = as.numeric(niter),
    B_min = floor(B_min),
    B_max = if (is.finite(B_max)) floor(B_max) else Inf,
    multiflip = isTRUE(multiflip),
    psingle = if (is.null(psingle)) -1 else as.numeric(psingle),
    psplit = as.numeric(psplit),
    pmerge = as.numeric(pmerge),
    pmergesplit = as.numeric(pmergesplit),
    gibbs_sweeps = as.numeric(max(1L, as.integer(gibbs_sweeps))),
    ls_shuffle = isTRUE(ls_shuffle)
  )
  if (is.null(c)) {
    opts$c <- if (nested) 1.0 else 0.5
  } else if (is.numeric(c) && length(c) == 1L && !is.na(c) && c >= 0) {
    opts$c <- as.numeric(c)
  } else if (nested && is.numeric(c) && length(c) >= 1L && all(!is.na(c) & c >= 0)) {
    opts$c <- as.numeric(c[1L])
    opts$c_levels <- as.numeric(c)
  } else {
    stop("c must be a non-negative number (or per-level vector for nested states)",
         call. = FALSE)
  }
  opts
}

run_native_mcmc <- function(state, mcmc_opts, eq_opts) {
  nested <- is_nested_state(state)
  weighted <- is_weighted_state(state)
  layered <- is_layered_state(state)
  pp <- is_pp_state(state)
  ranked <- is_ranked_state(state)
  bipartite <- is_bipartite_state(state)
  pointer <- graph_pointer(state$graph)
  dist <- identical(state$degree_prior, "distributed")
  if (bipartite) {
    res <- cpp_bipartite_mcmc(pointer, state$type, state$blocks,
                              state$degree_corrected, dist,
                              identical(state$edge_prior, "bipartite"),
                              state$weights, state$types, state$cpp_params,
                              mcmc_opts, eq_opts, q_cache_max())
    wspec <- list(W = state$weights, types = state$types,
                  params = state$params, cpp_params = state$cpp_params,
                  weighted = state$weighted)
    new_state <- new_bipartite_state(state$graph, state$type, res$blocks,
                                     state$degree_corrected, state$degree_prior,
                                     state$edge_prior, wspec, state$dim_labels,
                                     fit = list(description_length = res$entropy))
  } else if (pp) {
    res <- cpp_pp_mcmc(pointer, state$blocks, state$uniform,
                       state$degree_corrected, dist, mcmc_opts, eq_opts,
                       q_cache_max())
    new_state <- new_pp_state(state$graph, res$blocks, state$uniform,
                              state$degree_corrected, state$degree_prior,
                              fit = list(description_length = res$entropy))
  } else if (ranked && nested) {
    res <- cpp_nested_ranked_mcmc(pointer, state$levels, state$u,
                                  state$distance, state$unique,
                                  as.numeric(state$M), state$degree_corrected,
                                  dist, mcmc_opts, eq_opts, q_cache_max())
    new_state <- new_nested_ranked_state(state$graph, res$bs, res$u,
                                         state$distance, state$unique, state$M,
                                         state$degree_corrected,
                                         state$degree_prior,
                                         fit = list(description_length = res$entropy))
  } else if (ranked) {
    res <- cpp_ranked_mcmc(pointer, state$blocks, state$u,
                           state$distance, state$unique,
                           as.numeric(state$M), state$degree_corrected,
                           dist, mcmc_opts, eq_opts, q_cache_max())
    new_state <- new_ranked_state(state$graph, res$blocks, res$u,
                                  state$distance, state$unique, state$M,
                                  state$degree_corrected, state$degree_prior,
                                  fit = list(description_length = res$entropy))
  } else if (layered && nested) {
    w_mat <- if (isTRUE(state$weighted)) state$weights else matrix(0, nrow = 0L, ncol = 0L)
    w_types <- if (isTRUE(state$weighted)) state$types else character(0)
    w_params <- if (isTRUE(state$weighted)) state$cpp_params else list()
    res <- cpp_nested_layered_mcmc(pointer, state$levels, state$layers,
                                   state$independent, state$degree_corrected,
                                   dist, mcmc_opts, eq_opts, q_cache_max(),
                                   w_mat, w_types, w_params)
    linfo <- list(codes = state$layers, levels = state$layer_levels, C = state$C)
    wspec <- list(W = state$weights, types = state$types,
                  params = state$params, cpp_params = state$cpp_params,
                  weighted = isTRUE(state$weighted))
    new_state <- new_nested_layered_state(state$graph, res$bs, linfo,
                                          state$independent, state$degree_corrected,
                                          state$degree_prior,
                                          fit = list(description_length = res$entropy),
                                          wspec = wspec)
  } else if (layered) {
    w_mat <- if (isTRUE(state$weighted)) state$weights else matrix(0, nrow = 0L, ncol = 0L)
    w_types <- if (isTRUE(state$weighted)) state$types else character(0)
    w_params <- if (isTRUE(state$weighted)) state$cpp_params else list()
    res <- cpp_layered_mcmc(pointer, state$blocks, state$layers,
                            state$independent, state$degree_corrected,
                            dist, mcmc_opts, eq_opts, q_cache_max(),
                            w_mat, w_types, w_params)
    linfo <- list(codes = state$layers, levels = state$layer_levels, C = state$C)
    wspec <- list(W = state$weights, types = state$types,
                  params = state$params, cpp_params = state$cpp_params,
                  weighted = isTRUE(state$weighted))
    new_state <- new_layered_state(state$graph, res$blocks, linfo,
                                   state$independent, state$degree_corrected,
                                   state$degree_prior,
                                   fit = list(description_length = res$entropy),
                                   wspec = wspec)
  } else if (weighted && nested) {
    res <- cpp_nested_weighted_mcmc(pointer, state$levels, state$weights,
                                    state$types, state$cpp_params,
                                    state$degree_corrected, dist, mcmc_opts,
                                    eq_opts, q_cache_max())
    spec <- list(types = state$types, params = state$params,
                 cpp_params = state$cpp_params)
    new_state <- new_nested_weighted_state(state$graph, res$bs, state$weights,
                                           spec, state$degree_corrected,
                                           state$degree_prior,
                                           fit = list(description_length = res$entropy))
  } else if (weighted) {
    res <- cpp_weighted_mcmc(pointer, state$blocks, state$weights,
                             state$types, state$cpp_params,
                             state$degree_corrected, dist, mcmc_opts,
                             eq_opts, q_cache_max())
    spec <- list(types = state$types, params = state$params,
                 cpp_params = state$cpp_params)
    new_state <- new_weighted_state(state$graph, res$blocks, state$weights,
                                    spec, state$degree_corrected,
                                    state$degree_prior,
                                    fit = list(description_length = res$entropy))
  } else if (nested) {
    res <- cpp_nested_mcmc(pointer, state$levels, state$degree_corrected, dist,
                           mcmc_opts, eq_opts, q_cache_max())
    new_state <- new_nested_state(state$graph, res$bs, state$degree_corrected,
                                  state$degree_prior,
                                  fit = list(description_length = res$entropy))
  } else {
    res <- cpp_sbm_mcmc(pointer, match(state$blocks, state$labels),
                        state$degree_corrected, dist, mcmc_opts, eq_opts,
                        q_cache_max())
    new_state <- new_fitted_state(state$graph, res$blocks, state$degree_corrected,
                                  state$degree_prior, res$entropy)
  }
  out <- list(
    state = new_state,
    entropy = res$entropy,
    dS = res$dS,
    nattempts = res$nattempts,
    nmoves = res$nmoves,
    niter = res$niter
  )
  if (!is.null(res$history)) out$history <- res$history
  if (!is.null(res$samples)) out$samples <- res$samples
  if (!is.null(res$sample_entropies)) out$sample_entropies <- res$sample_entropies
  out
}

#' @noRd
#' @export
mcmc_sweep <- function(state, beta = 1, niter = 1L, c = NULL, d = 0.01,
                       multiflip = TRUE, B_min = 1L, B_max = Inf,
                       psingle = NULL, psplit = 1, pmerge = 1,
                       pmergesplit = 1, gibbs_sweeps = 10L,
                       ls_shuffle = TRUE) {
  validate_mcmc_state(state)
  mcmc_opts <- build_mcmc_opts(
    nested = is_nested_state(state), beta = beta, niter = niter, c = c, d = d,
    multiflip = multiflip, B_min = B_min, B_max = B_max, psingle = psingle,
    psplit = psplit, pmerge = pmerge, pmergesplit = pmergesplit,
    gibbs_sweeps = gibbs_sweeps, ls_shuffle = ls_shuffle
  )
  run_native_mcmc(state, mcmc_opts, list(mode = "sweep"))
}

#' @noRd
#' @export
mcmc_equilibrate <- function(state, wait = 1000L, nbreaks = 2L, max_niter = Inf,
                             force_niter = NULL, epsilon = 0, multiflip = TRUE,
                             mcmc_args = list(), history = FALSE,
                             callback = NULL, collect = FALSE, thin = 1L,
                             verbose = FALSE) {
  validate_mcmc_state(state)
  if (!is.list(mcmc_args)) stop("mcmc_args must be a named list", call. = FALSE)
  if (!is.numeric(wait) || length(wait) != 1L || is.na(wait) || wait < 1) {
    stop("wait must be a positive integer", call. = FALSE)
  }
  if (!is.numeric(nbreaks) || length(nbreaks) != 1L || is.na(nbreaks) || nbreaks < 1) {
    stop("nbreaks must be a positive integer", call. = FALSE)
  }
  if (!is.numeric(thin) || length(thin) != 1L || is.na(thin) || thin < 1) {
    stop("thin must be a positive integer", call. = FALSE)
  }
  mcmc_args$multiflip <- multiflip
  mcmc_opts <- do.call(build_mcmc_opts,
                       c(list(nested = is_nested_state(state)), mcmc_args))

  if (is.null(callback)) {
    eq_opts <- list(
      mode = "equilibrate",
      wait = floor(wait),
      nbreaks = floor(nbreaks),
      max_niter = if (is.finite(max_niter)) floor(max_niter) else Inf,
      force_niter = if (is.null(force_niter)) -1 else floor(force_niter),
      epsilon = as.numeric(epsilon),
      history = isTRUE(history),
      collect = isTRUE(collect),
      thin = floor(thin),
      verbose = isTRUE(verbose)
    )
    return(run_native_mcmc(state, mcmc_opts, eq_opts))
  }

  if (!is.function(callback)) stop("callback must be a function or NULL", call. = FALSE)
  count <- 0L
  break_count <- 0L
  niter <- 0L
  total_dS <- 0
  total_attempts <- 0
  total_moves <- 0
  cur <- state
  S <- if (is_nested_state(cur)) nested_description_length(cur) else block_description_length(cur)
  min_S <- max_S <- S
  m_eps <- 1e-6
  hist_rows <- list()
  samples <- list()
  cb_values <- list()

  while (count < wait) {
    step <- run_native_mcmc(cur, mcmc_opts, list(mode = "sweep"))
    cur <- step$state
    S <- step$entropy
    niter <- niter + 1L
    total_dS <- total_dS + step$dS
    total_attempts <- total_attempts + step$nattempts
    total_moves <- total_moves + step$nmoves

    if (!is.null(force_niter)) {
      max_S <- max(S, max_S)
      min_S <- min(S, min_S)
    } else {
      if (abs(step$dS) >= abs(S - step$dS) * epsilon) {
        if (S > max_S + m_eps) {
          max_S <- S
          count <- 0L
        } else if (S < min_S - m_eps) {
          min_S <- S
          count <- 0L
        } else {
          count <- count + 1L
        }
      } else {
        count <- count + 1L
      }
      if (count >= wait) {
        break_count <- break_count + 1L
        if (break_count < nbreaks) {
          count <- 0L
          min_S <- max_S <- S
        }
      }
    }

    extra <- callback(cur)
    if (!is.null(extra)) cb_values[[length(cb_values) + 1L]] <- extra
    if (isTRUE(collect) && (niter %% as.integer(thin) == 0L)) {
      samples[[length(samples) + 1L]] <- if (is_nested_state(cur)) {
        nested_levels(cur)
      } else {
        block_membership(cur)
      }
    }
    if (isTRUE(history)) {
      hist_rows[[niter]] <- data.frame(
        nattempts = step$nattempts, nmoves = step$nmoves,
        entropy = S, beta = mcmc_opts$beta
      )
    }
    if (!is.null(force_niter) && niter >= force_niter) break
    if (niter >= max_niter) break
  }

  out <- list(
    state = cur, entropy = S, dS = total_dS,
    nattempts = total_attempts, nmoves = total_moves, niter = niter
  )
  if (isTRUE(history)) out$history <- do.call(rbind, hist_rows)
  if (isTRUE(collect)) out$samples <- samples
  if (length(cb_values) > 0L) out$callback_values <- cb_values
  out
}

#' @noRd
#' @export
mcmc_anneal <- function(state, beta_range = c(1, 10), niter = 100L,
                        history = FALSE, mcmc_equilibrate_args = list(),
                        verbose = FALSE) {
  validate_mcmc_state(state)
  if (!is.numeric(beta_range) || length(beta_range) != 2L ||
      anyNA(beta_range) || any(beta_range <= 0) || any(!is.finite(beta_range))) {
    stop("beta_range must contain two finite positive numbers", call. = FALSE)
  }
  if (!is.numeric(niter) || length(niter) != 1L || is.na(niter) || niter < 1) {
    stop("niter must be a positive integer", call. = FALSE)
  }
  eq_args <- mcmc_equilibrate_args
  mcmc_args <- if (is.null(eq_args$mcmc_args)) list() else eq_args$mcmc_args
  mcmc_args$multiflip <- if (is.null(eq_args$multiflip)) TRUE else eq_args$multiflip
  mcmc_opts <- do.call(build_mcmc_opts,
                       c(list(nested = is_nested_state(state)), mcmc_args))
  eq_opts <- list(
    mode = "anneal",
    wait = if (is.null(eq_args$wait)) 100 else floor(eq_args$wait),
    nbreaks = if (is.null(eq_args$nbreaks)) 2 else floor(eq_args$nbreaks),
    max_niter = if (is.null(eq_args$max_niter)) Inf else eq_args$max_niter,
    force_niter = if (is.null(eq_args$force_niter)) -1 else floor(eq_args$force_niter),
    epsilon = if (is.null(eq_args$epsilon)) 0 else as.numeric(eq_args$epsilon),
    history = isTRUE(history),
    collect = FALSE,
    thin = 1,
    verbose = isTRUE(verbose),
    beta_start = as.numeric(beta_range[1L]),
    beta_end = as.numeric(beta_range[2L]),
    anneal_niter = floor(niter)
  )
  run_native_mcmc(state, mcmc_opts, eq_opts)
}

#' @noRd
#' @export
sample_posterior <- function(state, nsamples = 1000L, thin = 10L, burnin = NULL,
                             multiflip = TRUE, mcmc_args = list(),
                             equilibrate_args = list(), history = FALSE,
                             verbose = FALSE) {
  validate_mcmc_state(state)
  if (!is.numeric(nsamples) || length(nsamples) != 1L || is.na(nsamples) || nsamples < 1) {
    stop("nsamples must be a positive integer", call. = FALSE)
  }
  if (!is.numeric(thin) || length(thin) != 1L || is.na(thin) || thin < 1) {
    stop("thin must be a positive integer", call. = FALSE)
  }
  nsamples <- as.integer(nsamples)
  thin <- as.integer(thin)

  eq_res <- NULL
  cur <- state
  if (is.null(burnin)) {
    eq_call <- utils::modifyList(list(wait = 100L, nbreaks = 2L), equilibrate_args)
    eq_res <- do.call(mcmc_equilibrate, c(list(
      state = cur, multiflip = multiflip, mcmc_args = mcmc_args,
      history = history, verbose = verbose
    ), eq_call))
    cur <- eq_res$state
  } else if (is.numeric(burnin) && length(burnin) == 1L && burnin > 0) {
    eq_res <- mcmc_equilibrate(
      cur, force_niter = as.integer(burnin), multiflip = multiflip,
      mcmc_args = mcmc_args, history = history, verbose = verbose
    )
    cur <- eq_res$state
  }

  samp_res <- mcmc_equilibrate(
    cur, force_niter = nsamples * thin, multiflip = multiflip,
    mcmc_args = mcmc_args, history = history, collect = TRUE,
    thin = thin, verbose = verbose
  )
  list(
    state = samp_res$state,
    samples = samp_res$samples,
    entropies = samp_res$sample_entropies,
    nested = is_nested_state(state),
    equilibration = eq_res,
    sampling = samp_res[setdiff(names(samp_res), c("state", "samples", "sample_entropies"))]
  )
}
