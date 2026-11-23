validate_count <- function(n) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) ||
      !is.finite(n) || n < 0 || n > .Machine$integer.max || n != floor(n)) {
    stop("n must be a non-negative integer within R's integer range", call. = FALSE)
  }
  as.integer(n)
}

graph_pointer <- function(graph) {
  if (!inherits(graph, "graphtool_graph") || typeof(graph$pointer) != "externalptr") {
    stop("graph must be a graph created by graph_create()", call. = FALSE)
  }
  graph$pointer
}

graph_create <- function(n = 0L, directed = TRUE) {
  n <- validate_count(n)
  if (!is.logical(directed) || length(directed) != 1L || is.na(directed)) {
    stop("directed must be TRUE or FALSE", call. = FALSE)
  }
  structure(list(pointer = cpp_graph_create(n, directed), directed = directed),
            class = "graphtool_graph")
}

graph_add_vertices <- function(graph, n = 1L) {
  pointer <- graph_pointer(graph)
  n <- validate_count(n)
  if (cpp_num_vertices(pointer) + n > .Machine$integer.max) {
    stop("The graph cannot exceed R's integer vertex range", call. = FALSE)
  }
  cpp_add_vertices(pointer, n)
  invisible(graph)
}

graph_add_edges <- function(graph, edges, properties = NULL) {
  pointer <- graph_pointer(graph)
  count <- cpp_num_vertices(pointer)
  if (!is.matrix(edges) || ncol(edges) != 2L || !is.numeric(edges) ||
      anyNA(edges) || any(!is.finite(edges)) ||
      any(edges != floor(edges)) || any(edges < 1 | edges > count)) {
    stop("edges must be a two-column numeric matrix of valid 1-based vertex IDs",
         call. = FALSE)
  }
  if (!is.null(properties)) {
    if (is.data.frame(properties)) properties <- as.list(properties)
    if (!is.list(properties) || is.null(names(properties)) ||
        any(!nzchar(names(properties)))) {
      stop("properties must be a named list or data.frame of numeric edge vectors",
           call. = FALSE)
    }
    for (nm in names(properties)) {
      v <- properties[[nm]]
      if (!is.numeric(v) || length(v) != nrow(edges)) {
        stop(sprintf("edge property '%s' must be numeric with length nrow(edges)", nm),
             call. = FALSE)
      }
    }
  }
  storage.mode(edges) <- "integer"
  old_e <- as.integer(cpp_num_edges(pointer))
  cpp_add_edges(pointer, edges)
  if (!is.null(properties)) {
    for (nm in names(properties)) {
      cpp_set_edge_property(pointer, nm, as.numeric(properties[[nm]]),
                            TRUE, old_e)
    }
  }
  invisible(graph)
}

graph_set_edge_property <- function(graph, name, values) {
  pointer <- graph_pointer(graph)
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("name must be a non-empty character string", call. = FALSE)
  }
  if (!is.numeric(values) || length(values) != graph_num_edges(graph)) {
    stop("values must be a numeric vector of length graph_num_edges(graph)",
         call. = FALSE)
  }
  cpp_set_edge_property(pointer, name, as.numeric(values), FALSE, 0L)
  invisible(graph)
}

graph_edge_property <- function(graph, name) {
  pointer <- graph_pointer(graph)
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("name must be a non-empty character string", call. = FALSE)
  }
  cpp_get_edge_property(pointer, name)
}

graph_edge_properties <- function(graph) {
  cpp_edge_property_names(graph_pointer(graph))
}

transform_weights <- function(x, transform = c("identity", "log", "log1p",
                                               "sqrt", "sign", "abs",
                                               "rank", "zscore"),
                              shift = 0) {
  transform <- match.arg(transform)
  if (!is.numeric(x) || anyNA(x) || any(!is.finite(x))) {
    stop("x must be a finite numeric vector", call. = FALSE)
  }
  if (!is.numeric(shift) || length(shift) != 1L || is.na(shift) || !is.finite(shift)) {
    stop("shift must be a finite numeric scalar", call. = FALSE)
  }
  xs <- as.numeric(x) + shift
  log_jac <- NA_real_
  out <- switch(
    transform,
    identity = {
      log_jac <- 0
      xs
    },
    log = {
      if (any(xs <= 0)) stop("log transform requires strictly positive weights", call. = FALSE)
      log_jac <- -sum(log(xs))
      log(xs)
    },
    log1p = {
      if (any(xs <= -1)) stop("log1p transform requires weights > -1", call. = FALSE)
      log_jac <- -sum(log1p(xs))
      log1p(xs)
    },
    sqrt = {
      if (any(xs <= 0)) stop("sqrt transform requires strictly positive weights", call. = FALSE)
      log_jac <- sum(-log(2) - 0.5 * log(xs))
      sqrt(xs)
    },
    sign = {
      as.numeric(sign(xs))
    },
    abs = {
      abs(xs)
    },
    rank = {
      rank(xs, ties.method = "average")
    },
    zscore = {
      s <- stats::sd(xs)
      if (is.na(s) || s == 0) {
        stop("zscore transform requires nonzero standard deviation", call. = FALSE)
      }
      log_jac <- -length(xs) * log(s)
      (xs - mean(xs)) / s
    }
  )
  attr(out, "transform") <- transform
  attr(out, "log_jacobian") <- log_jac
  out
}

graph_num_vertices <- function(graph) {
  cpp_num_vertices(graph_pointer(graph))
}

graph_num_edges <- function(graph) {
  cpp_num_edges(graph_pointer(graph))
}

graph_edges <- function(graph) {
  cpp_edges(graph_pointer(graph))
}

graph_degree <- function(graph, mode = c("total", "in", "out")) {
  mode <- match.arg(mode)
  cpp_degree(graph_pointer(graph), mode)
}
