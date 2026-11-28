#' @noRd
#' @export
block_state <- function(graph, blocks = NULL, degree_corrected = TRUE) {
  pointer <- graph_pointer(graph)
  n <- cpp_num_vertices(pointer)
  if (!is.logical(degree_corrected) || length(degree_corrected) != 1L ||
      is.na(degree_corrected)) {
    stop("degree_corrected must be TRUE or FALSE", call. = FALSE)
  }
  if (is.null(blocks)) {
    blocks <- rep.int(1L, n)
  }
  if (!is.numeric(blocks) || length(blocks) != n || anyNA(blocks) ||
      any(!is.finite(blocks)) || any(blocks < 1) ||
      any(blocks > .Machine$integer.max) || any(blocks != floor(blocks))) {
    stop("blocks must contain one positive integer label per vertex", call. = FALSE)
  }
  blocks <- as.integer(blocks)
  labels <- sort(unique(blocks))
  index <- match(blocks, labels)
  statistics <- cpp_block_statistics(pointer, index, length(labels))
  statistics$num_edges <- cpp_num_edges(pointer)
  dimnames(statistics$counts) <- list(as.character(labels), as.character(labels))
  names(statistics$sizes) <- as.character(labels)
  structure(list(blocks = blocks, labels = labels, statistics = statistics,
                 graph = graph, degree_corrected = degree_corrected,
                 degree_prior = "distributed"),
            class = "rgraphtool_block_state")
}

validate_block_state <- function(state) {
  if (!inherits(state, "rgraphtool_block_state")) {
    stop("state must be created by block_state()", call. = FALSE)
  }
}

#' @noRd
#' @export
block_membership <- function(state) {
  if (inherits(state, "rgraphtool_weighted_block_state")) {
    validate_weighted_state(state)
    return(state$blocks)
  }
  if (inherits(state, "rgraphtool_layered_block_state")) {
    validate_layered_state(state)
    return(state$blocks)
  }
  if (inherits(state, "rgraphtool_pp_block_state")) {
    validate_pp_state(state)
    return(state$blocks)
  }
  if (inherits(state, "rgraphtool_ranked_block_state")) {
    validate_ranked_state(state)
    return(state$blocks)
  }
  if (inherits(state, "rgraphtool_bipartite_block_state")) {
    validate_bipartite_state(state)
    return(state$blocks)
  }
  if (inherits(state, "rgraphtool_modularity_state")) {
    return(state$blocks)
  }
  validate_block_state(state)
  state$blocks
}

#' @noRd
#' @export
block_sizes <- function(state) {
  if (inherits(state, "rgraphtool_weighted_block_state")) {
    validate_weighted_state(state)
    return(as.numeric(tabulate(state$blocks, nbins = max(state$blocks, 0L))))
  }
  if (inherits(state, "rgraphtool_layered_block_state")) {
    validate_layered_state(state)
    return(as.numeric(tabulate(state$blocks, nbins = max(state$blocks, 0L))))
  }
  if (inherits(state, "rgraphtool_pp_block_state")) {
    validate_pp_state(state)
    return(as.numeric(tabulate(state$blocks, nbins = max(state$blocks, 0L))))
  }
  if (inherits(state, "rgraphtool_ranked_block_state")) {
    validate_ranked_state(state)
    return(as.numeric(tabulate(state$blocks, nbins = max(state$blocks, 0L))))
  }
  if (inherits(state, "rgraphtool_bipartite_block_state")) {
    validate_bipartite_state(state)
    return(as.numeric(tabulate(state$blocks, nbins = max(state$blocks, 0L))))
  }
  if (inherits(state, "rgraphtool_modularity_state")) {
    return(as.numeric(tabulate(state$blocks, nbins = max(state$blocks, 0L))))
  }
  validate_block_state(state)
  state$statistics$sizes
}

#' @noRd
#' @export
block_edge_counts <- function(state) {
  validate_block_state(state)
  state$statistics$counts
}

#' @noRd
#' @export
block_log_likelihood <- function(state) {
  validate_block_state(state)
  stats <- state$statistics
  counts <- stats$counts
  sizes <- stats$sizes
  if (stats$directed) {
    numerator <- sum(lgamma(counts + 1))
    out <- rowSums(counts)
    incoming <- colSums(counts)
    if (state$degree_corrected) {
      denominator <- sum(lgamma(out + 1)) + sum(lgamma(incoming + 1))
    } else {
      denominator <- sum((out + incoming) * log(sizes))
    }
  } else {
    diagonal <- diag(counts)
    numerator <- sum(lgamma(counts[upper.tri(counts)] + 1)) +
      sum((diagonal / 2) * log(2) + lgamma(diagonal / 2 + 1))
    degree <- rowSums(counts)
    if (state$degree_corrected) {
      denominator <- sum(lgamma(degree + 1))
    } else {
      denominator <- sum(degree * log(sizes))
    }
  }
  numerator - denominator - stats$edge_factor +
    if (state$degree_corrected) stats$degree_factor else 0
}

block_partition_table <- function(state) {
  stats <- state$statistics
  max_degree <- if (stats$directed) sum(stats$counts) else sum(stats$counts) / 2
  if (!stats$directed) {
    max_degree <- max_degree * 2
  }
  n <- length(state$blocks)
  if (max_degree > .Machine$integer.max - 1 || n > .Machine$integer.max - 1) {
    stop("Graph is too large for the degree partition table", call. = FALSE)
  }
  cpp_log_partition_table(as.integer(max_degree), as.integer(n))
}

block_description_length_impl <- function(state, degree_prior, partition_table = NULL) {
  n <- length(state$blocks)
  if (n == 0L) {
    return(0)
  }
  stats <- state$statistics
  sizes <- stats$sizes
  counts <- stats$counts
  b <- length(sizes)
  edges <- sum(counts) / if (stats$directed) 1 else 2
  pairs <- if (stats$directed) b * b else b * (b + 1) / 2
  score <- -block_log_likelihood(state) +
    lchoose(n - 1, b - 1) + lgamma(n + 1) - sum(lgamma(sizes + 1)) +
    log(n) + lchoose(pairs + edges - 1, edges)
  if (!state$degree_corrected) {
    return(score)
  }

  out <- rowSums(counts)
  incoming <- colSums(counts)
  degree_sums <- if (stats$directed) list(out, incoming) else list(out)
  degree_vectors <- if (stats$directed) {
    list(stats$out_degrees, stats$in_degrees)
  } else {
    list(stats$out_degrees)
  }
  if (degree_prior == "uniform") {
    for (degrees in degree_sums) {
      score <- score + sum(lchoose(sizes + degrees - 1, degrees))
    }
  } else {
    if (is.null(partition_table)) {
      partition_table <- block_partition_table(state)
    }
    index <- match(state$blocks, state$labels)
    for (j in seq_along(degree_vectors)) {
      degrees <- degree_vectors[[j]]
      sums <- degree_sums[[j]]
      score <- score + sum(partition_table[cbind(sums + 1, sizes + 1)]) +
        sum(lgamma(sizes + 1))
      for (r in seq_along(sizes)) {
        frequencies <- table(degrees[index == r])
        score <- score - sum(lgamma(frequencies + 1))
      }
    }
  }
  score
}

#' @noRd
#' @export
block_description_length <- function(state, degree_prior = NULL) {
  validate_block_state(state)
  if (is.null(degree_prior)) degree_prior <- state$degree_prior
  degree_prior <- match.arg(degree_prior, c("distributed", "uniform"))
  # Graphs only grow, so an unchanged edge count means the snapshot is current.
  if (!is.null(state$graph) &&
      graph_num_edges(state$graph) == state$statistics$num_edges &&
      graph_num_vertices(state$graph) == length(state$blocks)) {
    return(cpp_sbm_entropy(graph_pointer(state$graph), match(state$blocks, state$labels),
                           state$degree_corrected, degree_prior == "distributed",
                           q_cache_max()))
  }
  block_description_length_impl(state, degree_prior)
}

q_cache_max <- function() {
  value <- getOption("rgraphtool.q_cache_max", 5000)
  if (!is.numeric(value) || length(value) != 1L || is.na(value) || value < 1) {
    stop("option rgraphtool.q_cache_max must be a positive number", call. = FALSE)
  }
  as.numeric(value)
}

multilevel_defaults <- list(
  r = 0.9, merge_sweeps = 10, mh_sweeps = 10, init_r = 0.99,
  init_min_iter = 8, init_min_iter_beta = 4, init_beta = 1, c = 0.5,
  epsilon = 1e-3, max_refine = 1000, bisection = TRUE, nthreads = 0L,
  verbose = FALSE
)

multilevel_control <- function(control, B_min = 1, B_max = Inf, nthreads = NULL) {
  if (!is.list(control) || (length(control) && is.null(names(control))) ||
      any(!names(control) %in% names(multilevel_defaults))) {
    stop("control must be a named list with entries among: ",
         paste(names(multilevel_defaults), collapse = ", "), call. = FALSE)
  }
  has_ctrl_thr <- !is.null(control$nthreads)
  control <- utils::modifyList(multilevel_defaults, control)
  for (name in setdiff(names(control), c("verbose", "bisection", "nthreads"))) {
    value <- control[[name]]
    if (!is.numeric(value) || length(value) != 1L || is.na(value) || value < 0) {
      stop("control$", name, " must be a non-negative number", call. = FALSE)
    }
  }
  if (!is.logical(control$bisection) || length(control$bisection) != 1L ||
      is.na(control$bisection)) {
    stop("control$bisection must be TRUE or FALSE", call. = FALSE)
  }
  if (control$r <= 0 || control$r >= 1) {
    stop("control$r must be strictly between 0 and 1", call. = FALSE)
  }
  if (!is.numeric(B_min) || length(B_min) != 1L || is.na(B_min) || B_min < 1 ||
      !is.numeric(B_max) || length(B_max) != 1L || is.na(B_max) || B_max < B_min) {
    stop("B_min and B_max must satisfy 1 <= B_min <= B_max", call. = FALSE)
  }
  control$B_min <- floor(B_min)
  control$B_max <- if (is.finite(B_max)) floor(B_max) else Inf
  control$nthreads <- if (!is.null(nthreads)) {
    resolve_nthreads(nthreads)
  } else if (has_ctrl_thr) {
    resolve_nthreads(control$nthreads)
  } else {
    resolve_nthreads(NULL)
  }
  control$verbose <- isTRUE(control$verbose)
  control
}

new_fitted_state <- function(graph, blocks, degree_corrected, degree_prior, entropy,
                             extra = list()) {
  state <- block_state(graph, blocks, degree_corrected)
  state$degree_prior <- degree_prior
  state$fit <- c(list(description_length = entropy), extra)
  state
}

#' @noRd
#' @export
minimize_blockmodel_dl <- function(graph, degree_corrected = TRUE,
                                   degree_prior = c("distributed", "uniform"),
                                   B_min = 1, B_max = Inf, refine = TRUE,
                                   epsilon = 1e-3, init_blocks = NULL,
                                   nthreads = NULL, control = list()) {
  pointer <- graph_pointer(graph)
  degree_prior <- match.arg(degree_prior)
  if (!is.logical(degree_corrected) || length(degree_corrected) != 1L ||
      is.na(degree_corrected)) {
    stop("degree_corrected must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(refine) || length(refine) != 1L || is.na(refine)) {
    stop("refine must be TRUE or FALSE", call. = FALSE)
  }
  n <- cpp_num_vertices(pointer)
  if (!is.null(init_blocks)) {
    if (inherits(init_blocks, "rgraphtool_block_state")) {
      validate_block_state(init_blocks)
      init_blocks <- init_blocks$blocks
    }
    init_blocks <- validate_label_vector(init_blocks, "init_blocks")
    if (length(init_blocks) != n) {
      stop("init_blocks must have length graph_num_vertices(graph)", call. = FALSE)
    }
  }
  control$epsilon <- epsilon
  control <- multilevel_control(control, B_min, B_max, nthreads = nthreads)
  result <- cpp_sbm_minimize(pointer, degree_corrected, degree_prior == "distributed",
                             control, refine, q_cache_max(), init_blocks)
  new_fitted_state(graph, result$blocks, degree_corrected, degree_prior, result$entropy,
                   list(multilevel_description_length = result$entropy_multilevel,
                        trace = result$trace))
}

#' @noRd
#' @export
refine_blockmodel_dl <- function(graph, state, degree_prior = NULL, epsilon = 1e-3,
                                 nthreads = NULL, control = list()) {
  pointer <- graph_pointer(graph)
  validate_block_state(state)
  if (length(state$blocks) != cpp_num_vertices(pointer)) {
    stop("graph and state must have the same number of vertices", call. = FALSE)
  }
  if (is.null(degree_prior)) degree_prior <- state$degree_prior
  degree_prior <- match.arg(degree_prior, c("distributed", "uniform"))
  control$epsilon <- epsilon
  control <- multilevel_control(control, nthreads = nthreads)
  result <- cpp_sbm_refine(pointer, match(state$blocks, state$labels),
                           state$degree_corrected, degree_prior == "distributed",
                           control, q_cache_max())
  new_fitted_state(graph, result$blocks, state$degree_corrected, degree_prior,
                   result$entropy)
}

#' @noRd
#' @export
compare_blockmodels <- function(graph, degree_prior = c("distributed", "uniform"),
                                nested = FALSE, nruns = 1L, cores = NULL, ...) {
  degree_prior <- match.arg(degree_prior)
  if (!is.logical(nested) || length(nested) != 1L || is.na(nested)) {
    stop("nested must be TRUE or FALSE", call. = FALSE)
  }
  nruns <- validate_count(nruns)
  if (nruns < 1L) stop("nruns must be positive", call. = FALSE)
  fit <- if (nested) minimize_nested_blockmodel_dl else minimize_blockmodel_dl
  score <- if (nested) nested_description_length else block_description_length
  states <- lapply(c(FALSE, TRUE), function(degree_corrected) {
    best <- NULL
    .starts <- parallel_starts(nruns, function(run) {
      state <- fit(graph, degree_corrected = degree_corrected,
                   degree_prior = degree_prior, ...)
      state
    }, cores = cores, graph = graph)
    for (run in seq_len(nruns)) {
      state <- .starts[[run]]
      if (is.null(best) || score(state) < score(best)) best <- state
    }
    best
  })
  names(states) <- c("non-degree-corrected", "degree-corrected")
  scores <- vapply(states, score, numeric(1))
  comparison <- data.frame(
    model = names(states),
    description_length = unname(scores),
    row.names = NULL
  )
  winner <- which.min(scores)
  structure(list(comparison = comparison,
                 log_posterior_odds = unname(scores[2] - scores[1]),
                 best_model = names(states)[winner],
                 best_state = states[[winner]],
                 states = states,
                 nested = nested),
            class = "rgraphtool_blockmodel_comparison")
}

#' @noRd
#' @export
print.rgraphtool_blockmodel_comparison <- function(x, ...) {
  cat(sprintf("<%s block model comparison>\n", if (isTRUE(x$nested)) "nested" else "flat"))
  for (i in seq_len(nrow(x$comparison))) {
    cat(sprintf("%-22s description length %12.4f nats\n",
                paste0(x$comparison$model[i], ":"), x$comparison$description_length[i]))
  }
  odds <- x$log_posterior_odds
  cat(sprintf("ln Lambda (NDC vs DC): %.4f, about 10^%.0f in favour of the %s model\n",
              odds, abs(odds) / log(10),
              if (odds > 0) "non-degree-corrected" else "degree-corrected"))
  invisible(x)
}
