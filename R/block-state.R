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
                 degree_corrected = degree_corrected,
                 degree_prior = "distributed"),
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

block_description_length <- function(state, degree_prior = NULL) {
  validate_block_state(state)
  if (is.null(degree_prior)) degree_prior <- state$degree_prior
  degree_prior <- match.arg(degree_prior, c("distributed", "uniform"))
  block_description_length_impl(state, degree_prior)
}

minimize_blockmodel_dl <- function(graph, degree_corrected = TRUE,
                                   degree_prior = c("distributed", "uniform"),
                                   restarts = 3L, max_sweeps = 10L) {
  degree_prior <- match.arg(degree_prior)
  restarts <- validate_count(restarts)
  max_sweeps <- validate_count(max_sweeps)
  if (restarts < 1L || max_sweeps < 1L) {
    stop("restarts and max_sweeps must be positive", call. = FALSE)
  }
  baseline <- block_state(graph, degree_corrected = degree_corrected)
  baseline$degree_prior <- degree_prior
  n <- length(baseline$blocks)
  if (n == 0L) {
    return(baseline)
  }
  partition_table <- if (degree_corrected && degree_prior == "distributed") {
    block_partition_table(baseline)
  } else {
    NULL
  }
  score_state <- function(blocks) {
    candidate <- block_state(graph, blocks, degree_corrected)
    candidate$degree_prior <- degree_prior
    list(state = candidate, score = block_description_length_impl(
      candidate, degree_prior, partition_table))
  }
  best <- list(state = baseline, score = block_description_length_impl(
    baseline, degree_prior, partition_table))
  for (restart in seq_len(restarts)) {
    groups <- max(1L, as.integer(ceiling(sqrt(n))))
    current <- if (restart == 1L) best else score_state(sample.int(groups, n, replace = TRUE))
    for (sweep in seq_len(max_sweeps)) {
      improved <- FALSE
      for (v in sample.int(n)) {
        blocks <- current$state$blocks
        candidates <- c(sort(unique(blocks)), max(blocks) + 1L)
        for (target in sample(candidates, length(candidates))) {
          if (target == blocks[v]) next
          proposal <- blocks
          proposal[v] <- target
          trial <- score_state(proposal)
          if (trial$score < current$score - 1e-10) {
            current <- trial
            blocks <- current$state$blocks
            improved <- TRUE
          }
        }
      }
      labels <- sort(unique(current$state$blocks))
      if (length(labels) > 1L) {
        for (i in seq_len(length(labels) - 1L)) {
          for (j in seq.int(i + 1L, length(labels))) {
            blocks <- current$state$blocks
            if (!any(blocks == labels[i]) || !any(blocks == labels[j])) next
            blocks[blocks == labels[j]] <- labels[i]
            trial <- score_state(blocks)
            if (trial$score < current$score - 1e-10) {
              current <- trial
              improved <- TRUE
            }
          }
        }
      }
      if (!improved) break
    }
    if (current$score < best$score) best <- current
  }
  best$state
}
