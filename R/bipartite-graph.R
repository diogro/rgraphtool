validate_bipartite_type <- function(graph, type = NULL) {
  graph_pointer(graph)
  nv <- as.integer(graph_num_vertices(graph))
  if (is.null(type)) {
    if ("type" %in% graph_vertex_properties(graph)) {
      type <- graph_vertex_property(graph, "type")
    } else {
      stop("type must be provided or stored as vertex property 'type' on graph",
           call. = FALSE)
    }
  }
  if (length(type) != nv || anyNA(type)) {
    stop("type must have length graph_num_vertices(graph) with no missing values",
         call. = FALSE)
  }
  if (is.logical(type)) {
    t_int <- as.integer(type)
  } else if (is.numeric(type)) {
    u <- sort(unique(type))
    if (length(u) > 2L || any(!is.finite(u))) {
      stop("bipartite type must take at most 2 distinct finite values", call. = FALSE)
    }
    if (all(u %in% c(0, 1))) {
      t_int <- as.integer(type)
    } else if (all(u %in% c(1, 2))) {
      t_int <- as.integer(type) - 1L
    } else {
      t_int <- match(type, u) - 1L
    }
  } else if (is.character(type) || is.factor(type)) {
    f <- as.factor(type)
    if (nlevels(f) > 2L) {
      stop("bipartite type must have at most 2 distinct levels", call. = FALSE)
    }
    t_int <- as.integer(f) - 1L
  } else {
    stop("type must be a logical, numeric, factor, or character vector", call. = FALSE)
  }
  as.integer(t_int)
}

get_dim_labels <- function(graph, dim_labels = NULL) {
  if (!is.null(dim_labels)) {
    if (!is.character(dim_labels) || length(dim_labels) != 2L ||
        anyNA(dim_labels) || any(!nzchar(dim_labels))) {
      stop("dim_labels must be a character vector of length 2", call. = FALSE)
    }
    return(unname(dim_labels))
  }
  if ("dim_labels" %in% graph_graph_properties(graph)) {
    dl <- graph_graph_property(graph, "dim_labels")
    if (is.character(dl) && length(dl) == 2L && !anyNA(dl) && all(nzchar(dl))) {
      return(unname(dl))
    }
  }
  c("row", "col")
}

bipartite_graph <- function(x, type = NULL, directed = FALSE,
                            dim_labels = NULL,
                            weight_mode = c("auto", "binary", "multigraph", "weight"),
                            weight_name = "weight",
                            covariates = NULL) {
  weight_mode <- match.arg(weight_mode)

  # Case 1: igraph object
  if (inherits(x, "igraph")) {
    g <- from_igraph(x)
    t_int <- validate_bipartite_type(g, type)
    check_bipartite_edges(g, t_int)
    graph_set_vertex_property(g, "type", as.logical(t_int))
    if (!is.null(dim_labels)) {
      graph_set_graph_property(g, "dim_labels", get_dim_labels(g, dim_labels))
    }
    return(g)
  }

  # Case 2: existing rgraphtool_graph
  if (inherits(x, "rgraphtool_graph")) {
    t_int <- validate_bipartite_type(x, type)
    check_bipartite_edges(x, t_int)
    graph_set_vertex_property(x, "type", as.logical(t_int))
    if (!is.null(dim_labels)) {
      graph_set_graph_property(x, "dim_labels", get_dim_labels(x, dim_labels))
    }
    return(x)
  }

  # Case 3: dense or sparse rectangular incidence matrix
  if (inherits(x, "Matrix")) {
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || (!is.numeric(x) && !is.logical(x))) {
    stop("x must be a rectangular numeric/logical matrix, rgraphtool_graph, or igraph object",
         call. = FALSE)
  }
  if (anyNA(x)) {
    stop("incidence matrix x must not contain NA values", call. = FALSE)
  }

  n1 <- nrow(x)
  n2 <- ncol(x)
  n <- n1 + n2

  if (is.null(dim_labels)) {
    dn <- names(dimnames(x))
    if (!is.null(dn) && length(dn) == 2L && all(!is.na(dn) & nzchar(dn))) {
      dim_labels <- unname(dn)
    } else {
      dim_labels <- c("row", "col")
    }
  } else {
    if (!is.character(dim_labels) || length(dim_labels) != 2L ||
        anyNA(dim_labels) || any(!nzchar(dim_labels))) {
      stop("dim_labels must be a character vector of length 2", call. = FALSE)
    }
    dim_labels <- unname(dim_labels)
  }

  vals_num <- if (is.logical(x)) {
    matrix(as.numeric(x), nrow = n1, ncol = n2)
  } else {
    x
  }

  if (weight_mode == "auto") {
    nz <- vals_num[vals_num != 0]
    if (!length(nz) || all(nz == 1)) {
      weight_mode <- "binary"
    } else if (all(nz > 0 & nz == floor(nz) & is.finite(nz))) {
      weight_mode <- "multigraph"
    } else {
      weight_mode <- "weight"
    }
  }

  g <- graph_create(n, directed = directed)
  graph_set_vertex_property(g, "type", c(rep(FALSE, n1), rep(TRUE, n2)))
  graph_set_graph_property(g, "dim_labels", dim_labels)

  rn <- rownames(x)
  cn <- colnames(x)
  if (!is.null(rn) || !is.null(cn)) {
    if (is.null(rn)) rn <- paste0(dim_labels[1L], "_", seq_len(n1))
    if (is.null(cn)) cn <- paste0(dim_labels[2L], "_", seq_len(n2))
    graph_set_vertex_property(g, "name", c(as.character(rn), as.character(cn)))
  }

  idx <- which(vals_num != 0, arr.ind = TRUE)
  if (nrow(idx) > 0L) {
    # Sort by row then col for deterministic edge ordering
    ord <- order(idx[, 1L], idx[, 2L])
    idx <- idx[ord, , drop = FALSE]
    wvals <- vals_num[idx]

    if (weight_mode == "binary") {
      edges <- cbind(as.integer(idx[, 1L]), as.integer(n1 + idx[, 2L]))
      props <- NULL
    } else if (weight_mode == "multigraph") {
      if (any(wvals < 0 | wvals != floor(wvals) | !is.finite(wvals))) {
        stop("weight_mode = 'multigraph' requires non-negative integer counts",
             call. = FALSE)
      }
      counts <- as.integer(wvals)
      edges <- cbind(
        rep.int(as.integer(idx[, 1L]), counts),
        rep.int(as.integer(n1 + idx[, 2L]), counts)
      )
      props <- NULL
    } else {
      edges <- cbind(as.integer(idx[, 1L]), as.integer(n1 + idx[, 2L]))
      props <- stats::setNames(list(as.numeric(wvals)), weight_name)
    }

    if (!is.null(covariates)) {
      if (!is.list(covariates) || is.null(names(covariates)) ||
          any(!nzchar(names(covariates)))) {
        stop("covariates must be a named list of matrices with dim(x)", call. = FALSE)
      }
      if (weight_mode == "multigraph") {
        stop("covariates cannot be combined with weight_mode = 'multigraph'; use weight_mode = 'weight' or 'binary'",
             call. = FALSE)
      }
      if (is.null(props)) props <- list()
      for (nm in names(covariates)) {
        cm <- covariates[[nm]]
        if (!is.matrix(cm) || !identical(dim(cm), c(n1, n2))) {
          stop(sprintf("covariate '%s' must be a %d x %d matrix", nm, n1, n2),
               call. = FALSE)
        }
        props[[nm]] <- as.numeric(cm[idx])
      }
    }

    graph_add_edges(g, edges, properties = props)
  }

  g
}

check_bipartite_edges <- function(graph, t_int) {
  el <- graph_edges(graph)
  if (nrow(el) > 0L) {
    bad <- which(t_int[el[, 1L]] == t_int[el[, 2L]])
    if (length(bad) > 0L) {
      stop(sprintf("Graph is not bipartite: found %d edge(s) connecting vertices of the same type (e.g., %d--%d)",
                   length(bad), el[bad[1L], 1L], el[bad[1L], 2L]), call. = FALSE)
    }
  }
  invisible(TRUE)
}

is_bipartite <- function(graph, type = NULL) {
  graph_pointer(graph)
  nv <- as.integer(graph_num_vertices(graph))
  if (nv == 0L) return(TRUE)

  if (!is.null(type) || "type" %in% graph_vertex_properties(graph)) {
    t_int <- tryCatch(validate_bipartite_type(graph, type), error = function(e) NULL)
    if (is.null(t_int)) return(FALSE)
    el <- graph_edges(graph)
    if (nrow(el) == 0L) return(TRUE)
    return(all(t_int[el[, 1L]] != t_int[el[, 2L]]))
  }

  # Fallback: BFS 2-coloring check
  el <- graph_edges(graph)
  if (nrow(el) == 0L) return(TRUE)
  if (any(el[, 1L] == el[, 2L])) return(FALSE)
  adj <- vector("list", nv)
  for (i in seq_len(nrow(el))) {
    u <- el[i, 1L]
    v <- el[i, 2L]
    adj[[u]] <- c(adj[[u]], v)
    adj[[v]] <- c(adj[[v]], u)
  }
  col <- rep.int(-1L, nv)
  for (s in seq_len(nv)) {
    if (col[s] != -1L) next
    col[s] <- 0L
    q <- s
    head <- 1L
    while (head <= length(q)) {
      u <- q[head]
      head <- head + 1L
      for (w in adj[[u]]) {
        if (col[w] == -1L) {
          col[w] <- 1L - col[u]
          q <- c(q, w)
        } else if (col[w] == col[u]) {
          return(FALSE)
        }
      }
    }
  }
  TRUE
}

as_incidence_matrix <- function(graph, type = NULL, weight = NULL) {
  t_int <- validate_bipartite_type(graph, type)
  check_bipartite_edges(graph, t_int)
  rows <- which(t_int == 0L)
  cols <- which(t_int == 1L)
  n1 <- length(rows)
  n2 <- length(cols)
  row_pos <- integer(length(t_int))
  col_pos <- integer(length(t_int))
  row_pos[rows] <- seq_len(n1)
  col_pos[cols] <- seq_len(n2)

  mat <- matrix(0, nrow = n1, ncol = n2)
  el <- graph_edges(graph)
  if (nrow(el) > 0L) {
    w <- if (is.null(weight)) rep(1, nrow(el)) else graph_edge_property(graph, weight)
    for (i in seq_len(nrow(el))) {
      u <- el[i, 1L]
      v <- el[i, 2L]
      if (t_int[u] == 0L) {
        mat[row_pos[u], col_pos[v]] <- mat[row_pos[u], col_pos[v]] + w[i]
      } else {
        mat[row_pos[v], col_pos[u]] <- mat[row_pos[v], col_pos[u]] + w[i]
      }
    }
  }

  if ("name" %in% graph_vertex_properties(graph)) {
    nms <- graph_vertex_property(graph, "name")
    if (is.character(nms) && length(nms) == length(t_int)) {
      rownames(mat) <- nms[rows]
      colnames(mat) <- nms[cols]
    }
  }
  dl <- get_dim_labels(graph)
  dimnames(mat) <- stats::setNames(
    list(rownames(mat), colnames(mat)),
    dl
  )
  mat
}
