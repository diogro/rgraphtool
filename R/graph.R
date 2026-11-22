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
  structure(list(pointer = cpp_graph_create(n, directed)), class = "graphtool_graph")
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

graph_add_edges <- function(graph, edges) {
  pointer <- graph_pointer(graph)
  count <- cpp_num_vertices(pointer)
  if (!is.matrix(edges) || ncol(edges) != 2L || !is.numeric(edges) ||
      anyNA(edges) || any(!is.finite(edges)) ||
      any(edges != floor(edges)) || any(edges < 1 | edges > count)) {
    stop("edges must be a two-column numeric matrix of valid 1-based vertex IDs",
         call. = FALSE)
  }
  storage.mode(edges) <- "integer"
  cpp_add_edges(pointer, edges)
  invisible(graph)
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
