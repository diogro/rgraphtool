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
  dimnames(statistics$counts) <- list(as.character(labels), as.character(labels))
  names(statistics$sizes) <- as.character(labels)
  structure(list(blocks = blocks, labels = labels, statistics = statistics,
                 degree_corrected = degree_corrected),
            class = "graphtool_block_state")
}

validate_block_state <- function(state) {
  if (!inherits(state, "graphtool_block_state")) {
    stop("state must be created by block_state()", call. = FALSE)
  }
}

block_membership <- function(state) {
  validate_block_state(state)
  state$blocks
}

block_sizes <- function(state) {
  validate_block_state(state)
  state$statistics$sizes
}

block_edge_counts <- function(state) {
  validate_block_state(state)
  state$statistics$counts
}

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
