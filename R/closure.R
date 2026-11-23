# Latent triadic closure blockmodel (LatentClosureBlockState, Peixoto 2022 PRX).
#
# Copyright (C) 2006-2026 Tiago de Paula Peixoto (graph-tool algorithms)
# Copyright (C) 2026 Diogo Melo (R port)
# LGPL-3+

latent_closure_block_state <- function(graph, L = 1L, nested = TRUE,
                                       blocks = NULL, gen = NULL,
                                       w0 = NULL, loops = NULL, ex = NULL,
                                       degree_corrected = TRUE,
                                       degree_prior = c("distributed", "uniform"),
                                       q_max = Inf) {
  graph_pointer(graph)
  if (isTRUE(graph$directed)) {
    stop("latent_closure_block_state() currently requires an undirected simple graph",
         call. = FALSE)
  }
  L <- as.integer(L)
  if (is.na(L) || L < 1L) stop("L must be an integer >= 1", call. = FALSE)
  degree_prior <- match.arg(degree_prior)
  dist_flag <- identical(degree_prior, "distributed")
  nested <- isTRUE(nested)
  degree_corrected <- isTRUE(degree_corrected)

  N <- graph_num_vertices(graph)
  raw_e <- graph_edges(graph)
  if (nrow(raw_e) > 0L) {
    if (any(raw_e[, 1L] == raw_e[, 2L])) {
      stop("latent_closure_block_state() expects a simple graph without self-loops", call. = FALSE)
    }
    u <- pmin(raw_e[, 1L], raw_e[, 2L])
    v <- pmax(raw_e[, 1L], raw_e[, 2L])
    obs_e <- cbind(as.integer(u), as.integer(v))
    if (anyDuplicated(paste0(u, ":", v))) {
      stop("latent_closure_block_state() expects a simple graph without parallel edges", call. = FALSE)
    }
  } else {
    obs_e <- matrix(integer(0), nrow = 0L, ncol = 2L)
  }
  m <- nrow(obs_e)

  gen_v <- if (is.null(gen)) rep(0L, m) else as.integer(gen)
  if (length(gen_v) != m || anyNA(gen_v) || any(gen_v < 0L | gen_v > L)) {
    stop(sprintf("gen must be an integer vector of length %d in 0..%d", m, L), call. = FALSE)
  }
  w0_v <- if (is.null(w0)) ifelse(gen_v == 0L, 1L, 0L) else as.integer(w0)
  if (length(w0_v) != m || anyNA(w0_v) || any(w0_v < 0L)) {
    stop(sprintf("w0 must be a non-negative integer vector of length %d", m), call. = FALSE)
  }
  loops_v <- if (is.null(loops)) rep(0L, N) else as.integer(loops)
  if (length(loops_v) != N || anyNA(loops_v) || any(loops_v < 0L)) {
    stop(sprintf("loops must be a non-negative integer vector of length %d", N), call. = FALSE)
  }
  if (is.null(ex)) {
    ex_list <- replicate(m, integer(0), simplify = FALSE)
  } else {
    if (!is.list(ex) || length(ex) != m) {
      stop(sprintf("ex must be a list of length %d", m), call. = FALSE)
    }
    ex_list <- lapply(ex, function(z) {
      if (length(z) == 0L) integer(0) else as.integer(z)
    })
  }

  bs <- .uncertain_validate_bs(N, blocks, nested)

  state <- structure(
    list(
      graph = graph,
      n_vertices = as.integer(N),
      L = L,
      nested = nested,
      degree_corrected = degree_corrected,
      degree_prior = degree_prior,
      dist = dist_flag,
      q_max = as.numeric(q_max),
      obs_edges = obs_e,
      gen = gen_v,
      w0 = w0_v,
      loops = loops_v,
      ex = ex_list,
      bs = bs
    ),
    class = "rgraphtool_latent_closure_block_state"
  )

  # Validate initial state immediately via C++
  invisible(closure_description_length(state, breakdown = FALSE))
  state
}

.check_closure_state <- function(state) {
  if (!inherits(state, "rgraphtool_latent_closure_block_state")) {
    stop("state must be created by latent_closure_block_state()", call. = FALSE)
  }
}

closure_description_length <- function(state, breakdown = FALSE) {
  .check_closure_state(state)
  res <- cpp_closure_entropy(
    state$n_vertices, state$L, state$obs_edges,
    state$gen, state$w0, state$loops, state$ex,
    state$bs, state$degree_corrected, state$dist, state$q_max
  )
  if (isTRUE(breakdown)) return(res)
  as.numeric(res$total)
}

.closure_init_sbm_bs <- function(state) {
  N <- state$n_vertices
  if (N <= 1L || nrow(state$obs_edges) == 0L) return(state$bs)
  if (length(unique(state$bs[[1L]])) > 1L) return(state$bs)
  g_sub <- graph_create(N, directed = FALSE)
  keep <- state$gen == 0L
  if (any(keep)) {
    graph_add_edges(g_sub, state$obs_edges[keep, , drop = FALSE])
  } else {
    return(state$bs)
  }
  if (state$nested) {
    fit <- minimize_nested_blockmodel_dl(
      g_sub, degree_corrected = state$degree_corrected,
      degree_prior = state$degree_prior
    )
    lapply(fit$levels, as.integer)
  } else {
    fit <- minimize_blockmodel_dl(
      g_sub, degree_corrected = state$degree_corrected,
      degree_prior = state$degree_prior
    )
    list(as.integer(fit$blocks))
  }
}

infer_triadic_closure <- function(state, nsamples = 200L, burnin = 100L, thin = 1L,
                                  beta = 1, pedges = 0.5, c = 1, d = 0.01,
                                  sample_loops = TRUE, init_sbm = TRUE,
                                  collect_partitions = FALSE) {
  .check_closure_state(state)
  nsamples <- as.integer(nsamples)
  burnin <- as.integer(burnin)
  thin <- as.integer(thin)
  if (is.na(nsamples) || nsamples < 1L) stop("nsamples must be >= 1", call. = FALSE)
  if (is.na(burnin) || burnin < 0L) stop("burnin must be >= 0", call. = FALSE)
  if (is.na(thin) || thin < 1L) stop("thin must be >= 1", call. = FALSE)

  bs_start <- if (isTRUE(init_sbm)) .closure_init_sbm_bs(state) else state$bs

  res <- cpp_closure_infer(
    state$n_vertices, state$L, state$obs_edges,
    state$gen, state$w0, state$loops, state$ex,
    bs_start, state$degree_corrected, state$dist,
    nsamples, burnin, thin, as.numeric(beta),
    as.numeric(pedges), as.numeric(c), as.numeric(d),
    isTRUE(sample_loops), isTRUE(collect_partitions), state$q_max
  )

  lp <- res$layer_probs
  colnames(lp) <- paste0("gen_", seq.int(0L, state$L))
  edges_df <- data.frame(
    from = state$obs_edges[, 1L],
    to = state$obs_edges[, 2L],
    substrate_prob = res$substrate_prob,
    closure_prob = res$closure_prob,
    as.data.frame(lp),
    stringsAsFactors = FALSE
  )

  vstats <- res$vertex_stats
  names(vstats) <- paste0("gen_", seq_len(state$L))

  new_state <- state
  new_state$gen <- as.integer(res$final_gen)
  new_state$w0 <- as.integer(res$final_w0)
  new_state$loops <- as.integer(res$final_loops)
  new_state$ex <- res$final_ex
  new_state$bs <- res$final_bs

  out <- list(
    edges = edges_df,
    vertex_stats = vstats,
    state = new_state,
    L = state$L,
    entropies = res$entropies,
    entropy = res$entropy,
    sbm_entropy = res$sbm_entropy,
    closure_entropy = res$closure_entropy,
    dS = res$dS,
    nattempts = res$nattempts,
    nmoves = res$nmoves
  )
  if (isTRUE(collect_partitions)) out$partitions <- res$partitions
  structure(out, class = "rgraphtool_closure_inference")
}

print.rgraphtool_latent_closure_block_state <- function(x, ...) {
  n_sub <- sum(x$gen == 0L)
  n_cl <- sum(x$gen > 0L)
  cat(sprintf("<rgraphtool_latent_closure_block_state (L=%d, %s)>\n",
              x$L, if (x$nested) "nested" else "flat"))
  cat(sprintf("  Vertices: %d | Edges: %d (substrate: %d, closure: %d)\n",
              x$n_vertices, nrow(x$obs_edges), n_sub, n_cl))
  cat(sprintf("  Blocks (level 1): %d | Degree-corrected: %s\n",
              length(unique(x$bs[[1L]])), if (x$degree_corrected) "TRUE" else "FALSE"))
  invisible(x)
}

print.rgraphtool_closure_inference <- function(x, ...) {
  cat(sprintf("<rgraphtool_closure_inference (L=%d)>\n", x$L))
  cat(sprintf("  Edges: %d | Mean substrate prob: %.3f | Mean closure prob: %.3f\n",
              nrow(x$edges), mean(x$edges$substrate_prob), mean(x$edges$closure_prob)))
  cat(sprintf("  Final DL: %.4f (SBM: %.4f, closure: %.4f) | Moves: %.0f / %.0f\n",
              x$entropy, x$sbm_entropy, x$closure_entropy, x$nmoves, x$nattempts))
  invisible(x)
}
