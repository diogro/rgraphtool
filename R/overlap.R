validate_half_edges <- function(graph, half_edges) {
  e <- graph_num_edges(graph)
  if (is.null(half_edges)) {
    out <- matrix(1L, nrow = e, ncol = 2L)
    colnames(out) <- c("source_block", "target_block")
    return(out)
  }
  half_edges <- as.matrix(half_edges)
  if (!is.numeric(half_edges) || nrow(half_edges) != e || ncol(half_edges) != 2L ||
      anyNA(half_edges) || any(!is.finite(half_edges)) || any(half_edges < 1) ||
      any(half_edges != floor(half_edges))) {
    stop("half_edges must be an E x 2 matrix of positive integer group labels",
         call. = FALSE)
  }
  storage.mode(half_edges) <- "integer"
  colnames(half_edges) <- c("source_block", "target_block")
  half_edges
}

new_overlap_state <- function(graph, half_edges, degree_corrected, degree_prior,
                              gen_overlap, overlap_info = NULL, fit = NULL) {
  if (is.null(overlap_info)) {
    overlap_info <- cpp_overlap_groups(graph_pointer(graph), half_edges)
  }
  structure(list(
    half_edges = half_edges,
    graph = graph,
    degree_corrected = degree_corrected,
    degree_prior = degree_prior,
    gen_overlap = gen_overlap,
    overlap = overlap_info,
    num_vertices = graph_num_vertices(graph),
    num_edges = graph_num_edges(graph),
    fit = fit
  ), class = "rgraphtool_overlap_block_state")
}

validate_overlap_state <- function(state) {
  if (!inherits(state, "rgraphtool_overlap_block_state")) {
    stop("state must be created by overlap_block_state() or minimize_overlap_blockmodel_dl()",
         call. = FALSE)
  }
  if (graph_num_edges(state$graph) != state$num_edges ||
      graph_num_vertices(state$graph) != state$num_vertices) {
    stop("the graph has changed since the overlap block state was created", call. = FALSE)
  }
}

#' @noRd
#' @export
overlap_block_state <- function(graph, half_edges = NULL, degree_corrected = TRUE,
                                degree_prior = c("distributed", "uniform"),
                                gen_overlap = TRUE) {
  graph_pointer(graph)
  degree_prior <- match.arg(degree_prior)
  if (!is.logical(degree_corrected) || length(degree_corrected) != 1L ||
      is.na(degree_corrected)) {
    stop("degree_corrected must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(gen_overlap) || length(gen_overlap) != 1L || is.na(gen_overlap)) {
    stop("gen_overlap must be TRUE or FALSE", call. = FALSE)
  }
  half_edges <- validate_half_edges(graph, half_edges)
  new_overlap_state(graph, half_edges, degree_corrected, degree_prior, gen_overlap)
}

#' @noRd
#' @export
overlap_description_length <- function(state, degree_prior = NULL) {
  validate_overlap_state(state)
  if (is.null(degree_prior)) degree_prior <- state$degree_prior
  degree_prior <- match.arg(degree_prior, c("distributed", "uniform"))
  cpp_overlap_entropy(graph_pointer(state$graph), state$half_edges,
                      state$degree_corrected, degree_prior == "distributed",
                      state$gen_overlap, q_cache_max())
}

#' @noRd
#' @export
overlap_groups <- function(state) {
  if (inherits(state, "rgraphtool_nested_overlap_block_state")) {
    validate_nested_overlap_state(state)
    return(state$overlap)
  }
  validate_overlap_state(state)
  state$overlap
}

#' @noRd
#' @export
overlap_majority_groups <- function(state) {
  overlap_groups(state)$majority
}

#' @noRd
#' @export
print.rgraphtool_overlap_block_state <- function(x, ...) {
  dl <- overlap_description_length(x)
  cat(sprintf("<overlap block state: %s, B = %d, description length %.4f nats>\n",
              if (x$degree_corrected) "degree-corrected" else "non-degree-corrected",
              x$overlap$B, dl))
  invisible(x)
}

#' @noRd
#' @export
minimize_overlap_blockmodel_dl <- function(graph, degree_corrected = TRUE,
                                           degree_prior = c("distributed", "uniform"),
                                           gen_overlap = TRUE, B_min = 1L,
                                           B_max = Inf, refine = TRUE,
                                           epsilon = 1e-3) {
  pointer <- graph_pointer(graph)
  degree_prior <- match.arg(degree_prior)
  if (!is.logical(degree_corrected) || length(degree_corrected) != 1L ||
      is.na(degree_corrected)) {
    stop("degree_corrected must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(gen_overlap) || length(gen_overlap) != 1L || is.na(gen_overlap)) {
    stop("gen_overlap must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.numeric(B_min) || length(B_min) != 1L || is.na(B_min) || B_min < 1 ||
      !is.numeric(B_max) || length(B_max) != 1L || is.na(B_max) || B_max < B_min) {
    stop("B_min and B_max must satisfy 1 <= B_min <= B_max", call. = FALSE)
  }
  # Seed with fast non-overlapping multilevel SBM fit as one of the starting points
  init_blocks <- rep.int(1L, graph_num_vertices(graph))
  if (graph_num_vertices(graph) > 1L && graph_num_edges(graph) > 0L) {
    flat_seed <- minimize_blockmodel_dl(graph, degree_corrected = degree_corrected,
                                        degree_prior = degree_prior,
                                        B_min = B_min, B_max = B_max,
                                        refine = TRUE)
    init_blocks <- match(block_membership(flat_seed), unique(block_membership(flat_seed)))
  }
  res <- cpp_overlap_minimize(pointer, as.integer(init_blocks), degree_corrected,
                              degree_prior == "distributed", gen_overlap,
                              as.integer(floor(B_min)),
                              if (is.finite(B_max)) floor(B_max) else Inf,
                              isTRUE(refine), as.numeric(epsilon), q_cache_max())
  new_overlap_state(graph, res$half_edges, degree_corrected, degree_prior,
                    gen_overlap, overlap_info = res$overlap,
                    fit = list(description_length = res$entropy))
}

# ---------------------------------------------------------------------------
# Nested overlapping block model (OverlapBlockState at level 0 + dense upper levels)

new_nested_overlap_state <- function(graph, half_edges, upper_levels,
                                     degree_corrected, degree_prior, gen_overlap,
                                     overlap_info = NULL, fit = NULL) {
  if (is.null(overlap_info)) {
    overlap_info <- cpp_overlap_groups(graph_pointer(graph), half_edges)
  }
  structure(list(
    half_edges = half_edges,
    upper_levels = upper_levels,
    graph = graph,
    degree_corrected = degree_corrected,
    degree_prior = degree_prior,
    gen_overlap = gen_overlap,
    overlap = overlap_info,
    num_vertices = graph_num_vertices(graph),
    num_edges = graph_num_edges(graph),
    fit = fit
  ), class = "rgraphtool_nested_overlap_block_state")
}

validate_nested_overlap_state <- function(state) {
  if (!inherits(state, "rgraphtool_nested_overlap_block_state")) {
    stop("state must be created by nested_overlap_block_state() or minimize_nested_overlap_blockmodel_dl()",
         call. = FALSE)
  }
  if (graph_num_edges(state$graph) != state$num_edges ||
      graph_num_vertices(state$graph) != state$num_vertices) {
    stop("the graph has changed since the nested overlap state was created", call. = FALSE)
  }
}

#' @noRd
#' @export
nested_overlap_block_state <- function(graph, half_edges, upper_levels,
                                       degree_corrected = TRUE,
                                       degree_prior = c("distributed", "uniform"),
                                       gen_overlap = TRUE) {
  graph_pointer(graph)
  degree_prior <- match.arg(degree_prior)
  half_edges <- validate_half_edges(graph, half_edges)
  info <- cpp_overlap_groups(graph_pointer(graph), half_edges)
  if (!is.list(upper_levels)) {
    stop("upper_levels must be a list of integer label vectors", call. = FALSE)
  }
  expected <- info$B
  clean_up <- vector("list", length(upper_levels))
  for (l in seq_along(upper_levels)) {
    b <- upper_levels[[l]]
    if (!is.numeric(b) || length(b) != expected || anyNA(b) || any(!is.finite(b)) ||
        any(b < 1) || any(b != floor(b))) {
      stop(sprintf("upper_levels[[%d]] must have %d positive integer labels", l, expected),
           call. = FALSE)
    }
    b <- match(b, unique(b))
    clean_up[[l]] <- as.integer(b)
    expected <- max(c(0L, b))
  }
  new_nested_overlap_state(graph, half_edges, clean_up, degree_corrected,
                           degree_prior, gen_overlap, overlap_info = info)
}

#' @noRd
#' @export
nested_overlap_level_entropy <- function(state, degree_prior = NULL) {
  validate_nested_overlap_state(state)
  if (is.null(degree_prior)) degree_prior <- state$degree_prior
  degree_prior <- match.arg(degree_prior, c("distributed", "uniform"))
  cpp_nested_overlap_entropy(graph_pointer(state$graph), state$half_edges,
                             state$upper_levels, state$degree_corrected,
                             degree_prior == "distributed", state$gen_overlap,
                             q_cache_max())
}

#' @noRd
#' @export
nested_overlap_description_length <- function(state, degree_prior = NULL) {
  sum(nested_overlap_level_entropy(state, degree_prior))
}

#' @noRd
#' @export
minimize_nested_overlap_blockmodel_dl <- function(graph, degree_corrected = TRUE,
                                                  degree_prior = c("distributed", "uniform"),
                                                  gen_overlap = TRUE, B_min = 1L,
                                                  B_max = Inf, refine = TRUE,
                                                  epsilon = 1e-3) {
  pointer <- graph_pointer(graph)
  degree_prior <- match.arg(degree_prior)
  if (graph_num_vertices(graph) == 0L) {
    stop("nested states need at least one vertex", call. = FALSE)
  }
  init_blocks <- rep.int(1L, graph_num_vertices(graph))
  if (graph_num_vertices(graph) > 1L && graph_num_edges(graph) > 0L) {
    nseed <- minimize_nested_blockmodel_dl(graph, degree_corrected = degree_corrected,
                                           degree_prior = degree_prior,
                                           B_min = B_min, B_max = B_max,
                                           refine = TRUE)
    init_blocks <- nested_levels(nseed)[[1L]]
  }
  res <- cpp_nested_overlap_minimize(pointer, as.integer(init_blocks),
                                     degree_corrected,
                                     degree_prior == "distributed", gen_overlap,
                                     as.integer(floor(B_min)),
                                     if (is.finite(B_max)) floor(B_max) else Inf,
                                     isTRUE(refine), as.numeric(epsilon), q_cache_max())
  new_nested_overlap_state(graph, res$half_edges, res$upper_levels,
                           degree_corrected, degree_prior, gen_overlap,
                           overlap_info = res$overlap,
                           fit = list(description_length = res$entropy))
}

#' @noRd
#' @export
print.rgraphtool_nested_overlap_block_state <- function(x, ...) {
  dl <- nested_overlap_description_length(x)
  cat(sprintf("<nested overlap block state: %s, %d levels, B_0 = %d, description length %.4f nats>\n",
              if (x$degree_corrected) "degree-corrected" else "non-degree-corrected",
              length(x$upper_levels) + 1L, x$overlap$B, dl))
  invisible(x)
}

#' @noRd
#' @export
overlap_mcmc_sweep <- function(state, beta = 1, niter = 1L, c = 0.5, d = 0.01,
                               thin = 1L, collect = FALSE, B_min = 1L,
                               B_max = Inf) {
  if (inherits(state, "rgraphtool_nested_overlap_block_state")) {
    validate_nested_overlap_state(state)
    res <- cpp_nested_overlap_mcmc(graph_pointer(state$graph), state$half_edges,
                                   state$upper_levels, state$degree_corrected,
                                   identical(state$degree_prior, "distributed"),
                                   state$gen_overlap, as.numeric(beta),
                                   as.numeric(c), as.numeric(d),
                                   as.integer(niter), as.integer(thin),
                                   isTRUE(collect), as.integer(floor(B_min)),
                                   if (is.finite(B_max)) floor(B_max) else Inf,
                                   q_cache_max())
    new_st <- new_nested_overlap_state(state$graph, res$half_edges,
                                       res$upper_levels, state$degree_corrected,
                                       state$degree_prior, state$gen_overlap,
                                       overlap_info = res$overlap,
                                       fit = list(description_length = res$entropy))
    out <- list(
      state = new_st,
      entropy = res$entropy,
      dS = res$dS,
      nattempts = res$nattempts,
      nmoves = res$nmoves
    )
    if (isTRUE(collect)) {
      out$majority_samples <- res$majority_samples
      out$entropies <- res$entropies
    }
    return(out)
  }
  validate_overlap_state(state)
  res <- cpp_overlap_mcmc(graph_pointer(state$graph), state$half_edges,
                          state$degree_corrected,
                          identical(state$degree_prior, "distributed"),
                          state$gen_overlap, as.numeric(beta), as.numeric(c),
                          as.numeric(d), as.integer(niter), as.integer(thin),
                          isTRUE(collect), as.integer(floor(B_min)),
                          if (is.finite(B_max)) floor(B_max) else Inf,
                          q_cache_max())
  new_st <- new_overlap_state(state$graph, res$half_edges, state$degree_corrected,
                              state$degree_prior, state$gen_overlap,
                              overlap_info = res$overlap,
                              fit = list(description_length = res$entropy))
  out <- list(
    state = new_st,
    entropy = res$entropy,
    dS = res$dS,
    nattempts = res$nattempts,
    nmoves = res$nmoves
  )
  if (isTRUE(collect)) {
    out$majority_samples <- res$majority_samples
    out$entropies <- res$entropies
  }
  out
}
