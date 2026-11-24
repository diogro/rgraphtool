validate_count <- function(n) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) ||
      !is.finite(n) || n < 0 || n > .Machine$integer.max || n != floor(n)) {
    stop("n must be a non-negative integer within R's integer range", call. = FALSE)
  }
  as.integer(n)
}

graph_pointer <- function(graph) {
  if (!inherits(graph, "rgraphtool_graph") || typeof(graph$pointer) != "externalptr") {
    stop("graph must be a graph created by graph_create()", call. = FALSE)
  }
  graph$pointer
}

graph_attrs_env <- function(graph) {
  graph_pointer(graph)
  env <- graph$attrs
  if (!is.environment(env)) {
    stop("graph has no attribute store; create graphs with graph_create()", call. = FALSE)
  }
  env
}

new_graph_attrs <- function() {
  env <- new.env(parent = emptyenv())
  env$graph_props <- list()
  env$vertex_props <- list()
  env$edge_props <- list()
  env
}

#' @noRd
#' @export
graph_create <- function(n = 0L, directed = TRUE) {
  n <- validate_count(n)
  if (!is.logical(directed) || length(directed) != 1L || is.na(directed)) {
    stop("directed must be TRUE or FALSE", call. = FALSE)
  }
  structure(list(pointer = cpp_graph_create(n, directed),
                 directed = directed,
                 attrs = new_graph_attrs()),
            class = "rgraphtool_graph")
}

#' @noRd
#' @export
graph_add_vertices <- function(graph, n = 1L) {
  pointer <- graph_pointer(graph)
  n <- validate_count(n)
  if (cpp_num_vertices(pointer) + n > .Machine$integer.max) {
    stop("The graph cannot exceed R's integer vertex range", call. = FALSE)
  }
  cpp_add_vertices(pointer, n)
  if (n > 0L && is.environment(graph$attrs) && length(graph$attrs$vertex_props)) {
    for (nm in names(graph$attrs$vertex_props)) {
      cur <- graph$attrs$vertex_props[[nm]]
      if (is.matrix(cur)) {
        pad <- matrix(NA, nrow = n, ncol = ncol(cur))
        storage.mode(pad) <- storage.mode(cur)
        graph$attrs$vertex_props[[nm]] <- rbind(cur, pad)
      } else if (is.list(cur)) {
        graph$attrs$vertex_props[[nm]] <- c(cur, vector("list", n))
      } else {
        pad <- rep(NA, n)
        storage.mode(pad) <- storage.mode(cur)
        graph$attrs$vertex_props[[nm]] <- c(cur, pad)
      }
    }
  }
  invisible(graph)
}

#' @noRd
#' @export
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
      stop("properties must be a named list or data.frame of edge vectors",
           call. = FALSE)
    }
    for (nm in names(properties)) {
      v <- properties[[nm]]
      vlen <- if (is.matrix(v)) nrow(v) else length(v)
      if (vlen != nrow(edges)) {
        stop(sprintf("edge property '%s' must have length nrow(edges)", nm),
             call. = FALSE)
      }
    }
  }
  storage.mode(edges) <- "integer"
  old_e <- as.integer(cpp_num_edges(pointer))
  old_ord <- if (old_e > 0L && is.environment(graph$attrs) && length(graph$attrs$edge_props)) {
    cpp_edge_order(pointer)
  } else {
    integer(0)
  }
  cpp_add_edges(pointer, edges)
  new_e <- as.integer(cpp_num_edges(pointer))
  new_ord <- if (new_e > 0L && (length(old_ord) > 0L || !is.null(properties))) {
    cpp_edge_order(pointer)
  } else {
    integer(0)
  }
  if (new_e > old_e && is.environment(graph$attrs) && length(graph$attrs$edge_props)) {
    added <- new_e - old_e
    for (nm in names(graph$attrs$edge_props)) {
      cur <- graph$attrs$edge_props[[nm]]
      # Convert existing values back to insertion order, append NAs, then map to new edge order
      if (is.matrix(cur)) {
        by_idx <- matrix(NA, nrow = new_e, ncol = ncol(cur))
        storage.mode(by_idx) <- storage.mode(cur)
        by_idx[old_ord, ] <- cur
        graph$attrs$edge_props[[nm]] <- by_idx[new_ord, , drop = FALSE]
      } else if (is.list(cur)) {
        by_idx <- vector("list", new_e)
        by_idx[old_ord] <- cur
        graph$attrs$edge_props[[nm]] <- by_idx[new_ord]
      } else {
        by_idx <- rep(NA, new_e)
        storage.mode(by_idx) <- storage.mode(cur)
        by_idx[old_ord] <- cur
        graph$attrs$edge_props[[nm]] <- by_idx[new_ord]
      }
    }
  }
  if (!is.null(properties)) {
    for (nm in names(properties)) {
      v <- properties[[nm]]
      if (is.numeric(v) && !is.matrix(v)) {
        cpp_set_edge_property(pointer, nm, as.numeric(v), TRUE, old_e)
      } else {
        env <- graph_attrs_env(graph)
        if (is.matrix(v)) {
          by_idx <- matrix(NA, nrow = new_e, ncol = ncol(v))
          storage.mode(by_idx) <- storage.mode(v)
          if (!is.null(env$edge_props[[nm]]) && old_e > 0L) {
            by_idx[new_ord, ] <- env$edge_props[[nm]]
          }
          by_idx[(old_e + 1L):new_e, ] <- v
          env$edge_props[[nm]] <- by_idx[new_ord, , drop = FALSE]
        } else if (is.list(v)) {
          by_idx <- vector("list", new_e)
          if (!is.null(env$edge_props[[nm]]) && old_e > 0L) {
            by_idx[new_ord] <- env$edge_props[[nm]]
          }
          by_idx[(old_e + 1L):new_e] <- v
          env$edge_props[[nm]] <- by_idx[new_ord]
        } else {
          by_idx <- rep(NA, new_e)
          storage.mode(by_idx) <- storage.mode(v)
          if (!is.null(env$edge_props[[nm]]) && old_e > 0L) {
            by_idx[new_ord] <- env$edge_props[[nm]]
          }
          by_idx[(old_e + 1L):new_e] <- v
          env$edge_props[[nm]] <- by_idx[new_ord]
        }
      }
    }
  }
  invisible(graph)
}

#' @noRd
#' @export
graph_set_edge_property <- function(graph, name, values) {
  pointer <- graph_pointer(graph)
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("name must be a non-empty character string", call. = FALSE)
  }
  ne <- graph_num_edges(graph)
  vlen <- if (is.matrix(values)) nrow(values) else length(values)
  if (vlen != ne) {
    stop("values must have length graph_num_edges(graph)", call. = FALSE)
  }
  if (is.numeric(values) && !is.matrix(values)) {
    cpp_set_edge_property(pointer, name, as.numeric(values), FALSE, 0L)
    if (is.environment(graph$attrs) && !is.null(graph$attrs$edge_props[[name]])) {
      graph$attrs$edge_props[[name]] <- NULL
    }
  } else if (is.logical(values) || is.character(values) || is.list(values) || is.matrix(values)) {
    env <- graph_attrs_env(graph)
    env$edge_props[[name]] <- values
  } else {
    stop("values must be a numeric, logical, character, matrix, or list edge property",
         call. = FALSE)
  }
  invisible(graph)
}

#' @noRd
#' @export
graph_edge_property <- function(graph, name) {
  pointer <- graph_pointer(graph)
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("name must be a non-empty character string", call. = FALSE)
  }
  if (name %in% cpp_edge_property_names(pointer)) {
    return(cpp_get_edge_property(pointer, name))
  }
  if (is.environment(graph$attrs) && name %in% names(graph$attrs$edge_props)) {
    return(graph$attrs$edge_props[[name]])
  }
  stop(sprintf("Edge property '%s' does not exist on graph", name), call. = FALSE)
}

#' @noRd
#' @export
graph_edge_properties <- function(graph) {
  cpp_names <- cpp_edge_property_names(graph_pointer(graph))
  extra_names <- if (is.environment(graph$attrs)) names(graph$attrs$edge_props) else character(0)
  unique(c(cpp_names, extra_names))
}

#' @noRd
#' @export
graph_set_vertex_property <- function(graph, name, values) {
  env <- graph_attrs_env(graph)
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("name must be a non-empty character string", call. = FALSE)
  }
  nv <- graph_num_vertices(graph)
  vlen <- if (is.matrix(values)) nrow(values) else length(values)
  if (vlen != nv) {
    stop("values must have length graph_num_vertices(graph)", call. = FALSE)
  }
  if (!(is.numeric(values) || is.logical(values) || is.character(values) ||
        is.list(values) || is.matrix(values))) {
    stop("values must be a numeric, logical, character, matrix, or list vertex property",
         call. = FALSE)
  }
  env$vertex_props[[name]] <- values
  invisible(graph)
}

#' @noRd
#' @export
graph_vertex_property <- function(graph, name) {
  env <- graph_attrs_env(graph)
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("name must be a non-empty character string", call. = FALSE)
  }
  if (!name %in% names(env$vertex_props)) {
    stop(sprintf("Vertex property '%s' does not exist on graph", name), call. = FALSE)
  }
  env$vertex_props[[name]]
}

#' @noRd
#' @export
graph_vertex_properties <- function(graph) {
  env <- graph_attrs_env(graph)
  nms <- names(env$vertex_props)
  if (is.null(nms)) character(0) else nms
}

#' @noRd
#' @export
graph_set_graph_property <- function(graph, name, value) {
  env <- graph_attrs_env(graph)
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("name must be a non-empty character string", call. = FALSE)
  }
  env$graph_props[[name]] <- value
  invisible(graph)
}

#' @noRd
#' @export
graph_graph_property <- function(graph, name) {
  env <- graph_attrs_env(graph)
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("name must be a non-empty character string", call. = FALSE)
  }
  if (!name %in% names(env$graph_props)) {
    stop(sprintf("Graph property '%s' does not exist on graph", name), call. = FALSE)
  }
  env$graph_props[[name]]
}

#' @noRd
#' @export
graph_graph_properties <- function(graph) {
  env <- graph_attrs_env(graph)
  nms <- names(env$graph_props)
  if (is.null(nms)) character(0) else nms
}

#' @noRd
#' @export
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

#' @noRd
#' @export
graph_num_vertices <- function(graph) {
  cpp_num_vertices(graph_pointer(graph))
}

#' @noRd
#' @export
graph_num_edges <- function(graph) {
  cpp_num_edges(graph_pointer(graph))
}

#' @noRd
#' @export
graph_edges <- function(graph) {
  cpp_edges(graph_pointer(graph))
}

#' @noRd
#' @export
graph_degree <- function(graph, mode = c("total", "in", "out")) {
  mode <- match.arg(mode)
  cpp_degree(graph_pointer(graph), mode)
}

format_prop_summary <- function(props) {
  if (!length(props)) return(NULL)
  tags <- vapply(names(props), function(nm) {
    v <- props[[nm]]
    tp <- if (is.matrix(v)) {
      "mat"
    } else if (inherits(v, "gt_python_object") || inherits(v, "gt_python_object_list")) {
      "pyobj"
    } else if (is.raw(v)) {
      "raw"
    } else if (is.list(v)) {
      "list"
    } else if (is.logical(v)) {
      "lgl"
    } else if (is.integer(v)) {
      "int"
    } else if (is.numeric(v)) {
      "dbl"
    } else if (is.character(v)) {
      "chr"
    } else {
      typeof(v)
    }
    sprintf("%s (%s)", nm, tp)
  }, character(1))
  paste(tags, collapse = ", ")
}

#' @noRd
#' @export
print.rgraphtool_graph <- function(x, max_edges = 20L, max_ascii = 12L,
                                   ascii = NULL, ...) {
  valid <- tryCatch({
    graph_num_vertices(x)
    TRUE
  }, error = function(e) FALSE)

  if (!valid) {
    cat("rgraphtool graph (invalid or deserialized external pointer)\n")
    return(invisible(x))
  }

  nv <- as.integer(graph_num_vertices(x))
  ne <- as.integer(graph_num_edges(x))
  dir_str <- if (isTRUE(x$directed)) "directed" else "undirected"
  cat(sprintf("rgraphtool graph (%s, %d %s, %d %s)\n",
              dir_str,
              nv, if (nv == 1L) "vertex" else "vertices",
              ne, if (ne == 1L) "edge" else "edges"))

  # Property map summaries
  if (is.environment(x$attrs)) {
    gp <- format_prop_summary(x$attrs$graph_props)
    if (!is.null(gp)) cat("  graph props: ", gp, "\n", sep = "")
    vp <- format_prop_summary(x$attrs$vertex_props)
    if (!is.null(vp)) cat("  vertex props: ", vp, "\n", sep = "")
  }
  ep_names <- graph_edge_properties(x)
  if (length(ep_names)) {
    ep_list <- stats::setNames(
      lapply(ep_names, function(nm) graph_edge_property(x, nm)),
      ep_names
    )
    cat("  edge props: ", format_prop_summary(ep_list), "\n", sep = "")
  }

  if (ne == 0L) {
    return(invisible(x))
  }

  edges <- graph_edges(x)
  # Use short vertex names if present
  vlabels <- as.character(seq_len(nv))
  if (is.environment(x$attrs) && "name" %in% names(x$attrs$vertex_props)) {
    cand <- x$attrs$vertex_props[["name"]]
    if (is.character(cand) && length(cand) == nv && !anyNA(cand) &&
        all(nzchar(cand)) && max(nchar(cand)) <= 12L) {
      vlabels <- cand
    }
  }

  # Compact ASCII adjacency matrix for small graphs
  show_ascii <- if (is.null(ascii)) (nv >= 1L && nv <= max_ascii) else isTRUE(ascii)
  if (show_ascii && nv >= 1L) {
    adj <- matrix(0L, nrow = nv, ncol = nv)
    for (i in seq_len(ne)) {
      u <- edges[i, 1L]
      v <- edges[i, 2L]
      adj[u, v] <- adj[u, v] + 1L
      if (!isTRUE(x$directed) && u != v) {
        adj[v, u] <- adj[v, u] + 1L
      }
    }
    lbl <- if (max(nchar(vlabels)) <= 3L) vlabels else as.character(seq_len(nv))
    w <- max(2L, max(nchar(lbl)))
    fmt_lbl <- paste0("%", w, "s")
    cat("  adjacency:\n")
    cat("    ", strrep(" ", w + 1L),
        paste(sprintf(fmt_lbl, lbl), collapse = " "), "\n", sep = "")
    for (i in seq_len(nv)) {
      cells <- ifelse(adj[i, ] == 0L, ".",
                      ifelse(adj[i, ] <= 9L, as.character(adj[i, ]), "*"))
      cat("    ", sprintf(fmt_lbl, lbl[i]), ":",
          paste(sprintf(fmt_lbl, cells), collapse = " "), "\n", sep = "")
    }
  }

  # Compact wrapped edge preview
  arrow <- if (isTRUE(x$directed)) "->" else "--"
  k <- min(ne, as.integer(max_edges))
  if (k > 0L) {
    parts <- sprintf("%s%s%s", vlabels[edges[seq_len(k), 1L]], arrow,
                     vlabels[edges[seq_len(k), 2L]])
    if (ne > k) {
      parts <- c(parts, sprintf("... (+%d more)", ne - k))
    }
    width <- max(40L, getOption("width", 80L) - 4L)
    lines <- strwrap(paste(parts, collapse = ", "), width = width)
    cat("  edges:\n")
    for (ln in lines) cat("    ", ln, "\n", sep = "")
  }

  invisible(x)
}
