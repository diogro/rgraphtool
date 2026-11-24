validate_hierarchy <- function(levels, n) {
  if (!is.list(levels) || length(levels) < 1L) {
    stop("levels must be a non-empty list of integer label vectors", call. = FALSE)
  }
  expected <- n
  out <- vector("list", length(levels))
  for (l in seq_along(levels)) {
    b <- levels[[l]]
    if (!is.numeric(b) || length(b) != expected || anyNA(b) || any(!is.finite(b)) ||
        any(b < 1) || any(b != floor(b))) {
      stop(sprintf(
        "levels[[%d]] must contain one positive integer label per %s (%d values)",
        l, if (l == 1L) "vertex" else sprintf("group of levels[[%d]]", l - 1L), expected
      ), call. = FALSE)
    }
    # relabel contiguously in order of first appearance
    b <- match(b, unique(b))
    out[[l]] <- as.integer(b)
    expected <- max(c(0L, b))
  }
  out
}

new_nested_state <- function(graph, levels, degree_corrected, degree_prior,
                             fit = NULL) {
  structure(list(levels = levels, graph = graph, degree_corrected = degree_corrected,
                 degree_prior = degree_prior, num_edges = graph_num_edges(graph),
                 fit = fit),
            class = "rgraphtool_nested_block_state")
}

validate_nested_state <- function(state) {
  if (!inherits(state, "rgraphtool_nested_block_state")) {
    stop("state must be created by nested_block_state() or minimize_nested_blockmodel_dl()",
         call. = FALSE)
  }
  if (graph_num_edges(state$graph) != state$num_edges ||
      graph_num_vertices(state$graph) != length(state$levels[[1]])) {
    stop("the graph has changed since the nested state was created", call. = FALSE)
  }
}

#' @noRd
#' @export
nested_block_state <- function(graph, levels, degree_corrected = TRUE,
                               degree_prior = c("distributed", "uniform")) {
  pointer <- graph_pointer(graph)
  degree_prior <- match.arg(degree_prior)
  if (!is.logical(degree_corrected) || length(degree_corrected) != 1L ||
      is.na(degree_corrected)) {
    stop("degree_corrected must be TRUE or FALSE", call. = FALSE)
  }
  n <- cpp_num_vertices(pointer)
  if (n == 0) {
    stop("nested states need at least one vertex", call. = FALSE)
  }
  new_nested_state(graph, validate_hierarchy(levels, n), degree_corrected, degree_prior)
}

#' @noRd
#' @export
nested_level_entropy <- function(state, degree_prior = NULL) {
  validate_nested_state(state)
  if (is.null(degree_prior)) degree_prior <- state$degree_prior
  degree_prior <- match.arg(degree_prior, c("distributed", "uniform"))
  cpp_nested_entropy(graph_pointer(state$graph), state$levels, state$degree_corrected,
                     degree_prior == "distributed", q_cache_max())
}

#' @noRd
#' @export
nested_description_length <- function(state, degree_prior = NULL) {
  sum(nested_level_entropy(state, degree_prior))
}

#' @noRd
#' @export
nested_levels <- function(state) {
  if (inherits(state, "rgraphtool_nested_weighted_block_state")) {
    validate_nested_weighted_state(state)
    return(state$levels)
  }
  if (inherits(state, "rgraphtool_nested_layered_block_state")) {
    validate_nested_layered_state(state)
    return(state$levels)
  }
  if (inherits(state, "rgraphtool_nested_ranked_block_state")) {
    validate_nested_ranked_state(state)
    return(state$levels)
  }
  validate_nested_state(state)
  state$levels
}

#' @noRd
#' @export
nested_project <- function(state, level = 0L) {
  if (inherits(state, "rgraphtool_nested_weighted_block_state")) {
    validate_nested_weighted_state(state)
  } else if (inherits(state, "rgraphtool_nested_layered_block_state")) {
    validate_nested_layered_state(state)
  } else if (inherits(state, "rgraphtool_nested_ranked_block_state")) {
    validate_nested_ranked_state(state)
  } else {
    validate_nested_state(state)
  }
  depth <- length(state$levels)
  if (!is.numeric(level) || length(level) != 1L || is.na(level) || level < 0 ||
      level >= depth || level != floor(level)) {
    stop(sprintf("level must be an integer between 0 and %d", depth - 1L), call. = FALSE)
  }
  b <- state$levels[[1]]
  for (l in seq_len(level)) {
    b <- state$levels[[l + 1L]][b]
  }
  b
}

#' @noRd
#' @export
nested_level_state <- function(state, level = 0L) {
  b <- nested_project(state, level)
  out <- block_state(state$graph, b, state$degree_corrected)
  out$degree_prior <- state$degree_prior
  out
}

#' @noRd
#' @export
nested_summary <- function(state) {
  if (inherits(state, "rgraphtool_nested_weighted_block_state")) {
    validate_nested_weighted_state(state)
    lv_ent <- nested_weighted_level_entropy(state)
  } else if (inherits(state, "rgraphtool_nested_layered_block_state")) {
    validate_nested_layered_state(state)
    lv_ent <- nested_layered_level_entropy(state)
  } else if (inherits(state, "rgraphtool_nested_ranked_block_state")) {
    validate_nested_ranked_state(state)
    lv_ent <- nested_ranked_level_entropy(state)
  } else {
    validate_nested_state(state)
    lv_ent <- nested_level_entropy(state)
  }
  sizes <- c(length(state$levels[[1]]), vapply(state$levels, max, integer(1)))
  depth <- length(state$levels)
  out <- data.frame(level = seq_len(depth) - 1L, N = sizes[seq_len(depth)],
                    B = sizes[-1L], description_length = lv_ent)
  last <- match(1L, out$N)
  if (!is.na(last)) out <- out[seq_len(last), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' @noRd
#' @export
print.rgraphtool_nested_block_state <- function(x, ...) {
  summary <- nested_summary(x)
  cat(sprintf("<nested block state: %s, %d levels, description length %.4f nats>\n",
              if (x$degree_corrected) "degree-corrected" else "non-degree-corrected",
              length(x$levels), sum(nested_level_entropy(x))))
  for (i in seq_len(nrow(summary))) {
    cat(sprintf("l: %d, N: %d, B: %d\n", summary$level[i], summary$N[i], summary$B[i]))
  }
  invisible(x)
}

#' @noRd
#' @export
minimize_nested_blockmodel_dl <- function(graph, degree_corrected = TRUE,
                                          degree_prior = c("distributed", "uniform"),
                                          B_min = 1, B_max = Inf, refine = TRUE,
                                          epsilon = 1e-3, simple_init = FALSE,
                                          simple_init_reverse = FALSE,
                                          base_factor = 10, B_min_base = 1000L,
                                          top_factor = 2, strict = FALSE,
                                          init_levels = NULL, nthreads = NULL,
                                          control = list()) {
  pointer <- graph_pointer(graph)
  degree_prior <- match.arg(degree_prior)
  if (!is.logical(degree_corrected) || length(degree_corrected) != 1L ||
      is.na(degree_corrected)) {
    stop("degree_corrected must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(refine) || length(refine) != 1L || is.na(refine)) {
    stop("refine must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(simple_init) || length(simple_init) != 1L || is.na(simple_init)) {
    stop("simple_init must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(simple_init_reverse) || length(simple_init_reverse) != 1L ||
      is.na(simple_init_reverse)) {
    stop("simple_init_reverse must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("strict must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.numeric(base_factor) || length(base_factor) != 1L || is.na(base_factor) ||
      !is.finite(base_factor) || base_factor <= 0) {
    stop("base_factor must be a positive finite number", call. = FALSE)
  }
  if (!is.numeric(B_min_base) || length(B_min_base) != 1L || is.na(B_min_base) ||
      !is.finite(B_min_base) || B_min_base < 1 || B_min_base != floor(B_min_base)) {
    stop("B_min_base must be a positive integer", call. = FALSE)
  }
  if (!is.numeric(top_factor) || length(top_factor) != 1L || is.na(top_factor) ||
      !is.finite(top_factor) || top_factor <= 1) {
    stop("top_factor must be a finite number greater than 1", call. = FALSE)
  }
  n <- cpp_num_vertices(pointer)
  if (n == 0) {
    stop("nested states need at least one vertex", call. = FALSE)
  }
  if (!is.null(init_levels)) {
    if (inherits(init_levels, "rgraphtool_nested_block_state")) {
      validate_nested_state(init_levels)
      init_levels <- init_levels$levels
    }
    init_levels <- validate_hierarchy(init_levels, n)
  }
  if (is.null(control$c)) control$c <- 1
  control$epsilon <- epsilon
  control <- multilevel_control(control, B_min, B_max, nthreads = nthreads)
  result <- cpp_nested_minimize(pointer, degree_corrected, degree_prior == "distributed",
                                control, refine, q_cache_max(),
                                isTRUE(simple_init), isTRUE(simple_init_reverse),
                                as.numeric(base_factor), as.integer(B_min_base),
                                as.numeric(top_factor), isTRUE(strict),
                                init_levels)
  new_nested_state(graph, result$bs, degree_corrected, degree_prior,
                   fit = list(description_length = sum(result$level_entropy),
                              trace = result$trace))
}
