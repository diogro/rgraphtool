validate_label_vector <- function(x, name = "x") {
  if (!is.numeric(x) || any(!is.na(x) & (!is.finite(x) | x < 1 | x != floor(x)))) {
    stop(name, " must contain positive integer group labels (or NA for absent nodes)",
         call. = FALSE)
  }
  as.integer(x)
}

validate_nested_label_list <- function(x, name = "x") {
  if (!is.list(x) || length(x) < 1L) {
    stop(name, " must be a non-empty list of integer label vectors", call. = FALSE)
  }
  lapply(seq_along(x), function(l) validate_label_vector(x[[l]], sprintf("%s[[%d]]", name, l)))
}

infer_nested_flag <- function(partitions, nested = NULL) {
  if (!is.list(partitions) || length(partitions) < 1L) {
    stop("partitions must be a non-empty list of partitions", call. = FALSE)
  }
  if (is.null(nested)) {
    return(is.list(partitions[[1L]]))
  }
  if (!is.logical(nested) || length(nested) != 1L || is.na(nested)) {
    stop("nested must be TRUE, FALSE, or NULL", call. = FALSE)
  }
  nested
}

normalize_partition_list <- function(partitions, nested) {
  if (nested) {
    lapply(seq_along(partitions), function(i) {
      validate_nested_label_list(partitions[[i]], sprintf("partitions[[%d]]", i))
    })
  } else {
    n <- length(partitions[[1L]])
    lapply(seq_along(partitions), function(i) {
      v <- validate_label_vector(partitions[[i]], sprintf("partitions[[%d]]", i))
      if (length(v) != n) {
        stop("all flat partitions must have the same length", call. = FALSE)
      }
      v
    })
  }
}

align_partition_labels <- function(x, y) {
  x <- validate_label_vector(x, "x")
  y <- validate_label_vector(y, "y")
  cpp_align_partition_labels(x, y)
}

partition_overlap <- function(x, y, norm = TRUE) {
  x <- validate_label_vector(x, "x")
  y <- validate_label_vector(y, "y")
  if (!is.logical(norm) || length(norm) != 1L || is.na(norm)) {
    stop("norm must be TRUE or FALSE", call. = FALSE)
  }
  cpp_partition_overlap(x, y, norm)
}

align_nested_partition_labels <- function(x, y) {
  x <- validate_nested_label_list(x, "x")
  y <- validate_nested_label_list(y, "y")
  cpp_align_nested_partition_labels(x, y)
}

nested_partition_overlap <- function(x, y, norm = TRUE) {
  x <- validate_nested_label_list(x, "x")
  y <- validate_nested_label_list(y, "y")
  if (!is.logical(norm) || length(norm) != 1L || is.na(norm)) {
    stop("norm must be TRUE or FALSE", call. = FALSE)
  }
  cpp_nested_partition_overlap(x, y, norm)
}

wrap_partition_mode <- function(raw) {
  structure(raw, class = "graphtool_partition_mode")
}

partition_mode_state <- function(partitions, relabel = TRUE, nested = NULL,
                                 converge = TRUE) {
  if (is.list(partitions) && !is.null(partitions$samples) &&
      !is.null(partitions$nested)) {
    if (is.null(nested)) nested <- partitions$nested
    partitions <- partitions$samples
  }
  nested <- infer_nested_flag(partitions, nested)
  partitions <- normalize_partition_list(partitions, nested)
  res <- cpp_partition_mode(partitions, isTRUE(relabel), nested, isTRUE(converge))
  wrap_partition_mode(res)
}

validate_partition_mode <- function(mode) {
  if (!inherits(mode, "graphtool_partition_mode")) {
    stop("mode must be created by partition_mode_state()", call. = FALSE)
  }
}

mode_marginals <- function(mode, probability = FALSE, level = 0L) {
  validate_partition_mode(mode)
  if (mode$nested) {
    if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
        level < 0 || level >= length(mode$marginals) || level != floor(level)) {
      stop("level out of range", call. = FALSE)
    }
    mat <- mode$marginals[[as.integer(level) + 1L]]
  } else {
    mat <- mode$marginals
  }
  if (isTRUE(probability)) {
    rs <- rowSums(mat)
    probs <- mat / pmax(rs, 1)
    probs[rs == 0, ] <- 0
    return(probs)
  }
  mat
}

mode_consensus <- function(mode) {
  validate_partition_mode(mode)
  mode$consensus
}

mode_uncertainty <- function(mode, mle = TRUE) {
  validate_partition_mode(mode)
  if (isTRUE(mle)) mode$posterior_cdev else mode$posterior_cdev_mean
}

mode_entropy <- function(mode, mle = TRUE) {
  validate_partition_mode(mode)
  if (isTRUE(mle)) mode$posterior_entropy else mode$posterior_entropy_mean
}

print.graphtool_partition_mode <- function(x, ...) {
  if (isTRUE(x$nested)) {
    cat(sprintf("<nested partition mode: M = %d, B = %s, uncertainty = %.4f, DL = %.4f nats>\n",
                x$M, paste(x$B, collapse = "-"), x$posterior_cdev, x$entropy))
  } else {
    cat(sprintf("<partition mode: M = %d, B = %d, uncertainty = %.4f, DL = %.4f nats>\n",
                x$M, x$B, x$posterior_cdev, x$entropy))
  }
  invisible(x)
}

mode_cluster_state <- function(partitions, nested = NULL, b = NULL, relabel = TRUE,
                               niter = 20L, converge = TRUE) {
  if (is.list(partitions) && !is.null(partitions$samples) &&
      !is.null(partitions$nested)) {
    if (is.null(nested)) nested <- partitions$nested
    partitions <- partitions$samples
  }
  nested <- infer_nested_flag(partitions, nested)
  partitions <- normalize_partition_list(partitions, nested)
  if (!is.null(b)) {
    b <- validate_label_vector(b, "b")
    if (length(b) != length(partitions)) {
      stop("b must have one cluster label per partition", call. = FALSE)
    }
  }
  if (!is.numeric(niter) || length(niter) != 1L || is.na(niter) || niter < 0 ||
      niter != floor(niter)) {
    stop("niter must be a non-negative integer", call. = FALSE)
  }
  res <- cpp_mode_cluster(partitions, nested, isTRUE(relabel), b,
                          as.integer(niter), isTRUE(converge))
  res$modes <- lapply(res$modes, wrap_partition_mode)
  res$nested <- nested
  structure(res, class = "graphtool_mode_cluster")
}

mode_summary <- function(cluster) {
  if (!inherits(cluster, "graphtool_mode_cluster")) {
    stop("cluster must be created by mode_cluster_state()", call. = FALSE)
  }
  modes <- cluster$modes
  data.frame(
    mode = seq_along(modes),
    M = vapply(modes, function(m) m$M, integer(1)),
    weight = as.numeric(cluster$weights),
    B = vapply(modes, function(m) if (isTRUE(m$nested)) m$B[1L] else m$B, integer(1)),
    uncertainty = vapply(modes, function(m) m$posterior_cdev, numeric(1)),
    posterior_entropy = vapply(modes, function(m) m$posterior_entropy, numeric(1)),
    description_length = vapply(modes, function(m) m$entropy, numeric(1)),
    row.names = NULL
  )
}

print.graphtool_mode_cluster <- function(x, ...) {
  cat(sprintf("<%s mode cluster state: K = %d modes (Be = %.2f), M = %d partitions, DL = %.4f nats>\n",
              if (isTRUE(x$nested)) "nested" else "flat",
              x$K, x$Be, length(x$clusters), x$entropy))
  print(mode_summary(x), row.names = FALSE)
  invisible(x)
}
