resolve_edge_weights_opt <- function(graph, weights = NULL) {
  graph_pointer(graph)
  if (is.null(weights)) return(NULL)
  if (is.character(weights) && length(weights) == 1L && !is.na(weights) && nzchar(weights)) {
    return(as.numeric(graph_edge_property(graph, weights)))
  }
  if (is.numeric(weights)) {
    if (length(weights) != graph_num_edges(graph)) {
      stop("weights must have length graph_num_edges(graph)", call. = FALSE)
    }
    return(as.numeric(weights))
  }
  stop("weights must be NULL, an edge property name, or a numeric vector of length E",
       call. = FALSE)
}

# ===========================================================================
# 1. Centrality Measures
# ===========================================================================

pagerank <- function(graph, damping = 0.85, pers = NULL, weights = NULL,
                     tol = 1e-6, max_iter = 1000L) {
  pointer <- graph_pointer(graph)
  if (!is.numeric(damping) || length(damping) != 1L || is.na(damping) ||
      damping < 0 || damping > 1) {
    stop("damping must be a number in [0, 1]", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol <= 0) {
    stop("tol must be a positive number", call. = FALSE)
  }
  max_iter <- validate_count(max_iter)
  w <- resolve_edge_weights_opt(graph, weights)
  res <- cpp_pagerank(pointer, as.numeric(damping), pers, w, as.numeric(tol), max_iter)
  out <- res$vector
  attr(out, "iterations") <- res$iterations
  out
}

betweenness <- function(graph, weights = NULL, normalized = TRUE) {
  pointer <- graph_pointer(graph)
  if (!is.logical(normalized) || length(normalized) != 1L || is.na(normalized)) {
    stop("normalized must be TRUE or FALSE", call. = FALSE)
  }
  w <- resolve_edge_weights_opt(graph, weights)
  cpp_betweenness(pointer, w, normalized)
}

closeness <- function(graph, weights = NULL, normalized = TRUE, harmonic = FALSE) {
  pointer <- graph_pointer(graph)
  if (!is.logical(normalized) || length(normalized) != 1L || is.na(normalized)) {
    stop("normalized must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(harmonic) || length(harmonic) != 1L || is.na(harmonic)) {
    stop("harmonic must be TRUE or FALSE", call. = FALSE)
  }
  w <- resolve_edge_weights_opt(graph, weights)
  cpp_closeness(pointer, w, normalized, harmonic)
}

eigenvector_centrality <- function(graph, weights = NULL, tol = 1e-6, max_iter = 1000L) {
  pointer <- graph_pointer(graph)
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol <= 0) {
    stop("tol must be a positive number", call. = FALSE)
  }
  max_iter <- validate_count(max_iter)
  w <- resolve_edge_weights_opt(graph, weights)
  cpp_eigenvector(pointer, w, as.numeric(tol), max_iter)
}

katz_centrality <- function(graph, alpha = 0.01, beta = NULL, weights = NULL,
                            normalized = TRUE, tol = 1e-6, max_iter = 1000L) {
  pointer <- graph_pointer(graph)
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || !is.finite(alpha)) {
    stop("alpha must be a finite scalar", call. = FALSE)
  }
  if (!is.logical(normalized) || length(normalized) != 1L || is.na(normalized)) {
    stop("normalized must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol <= 0) {
    stop("tol must be a positive number", call. = FALSE)
  }
  max_iter <- validate_count(max_iter)
  w <- resolve_edge_weights_opt(graph, weights)
  cpp_katz(pointer, as.numeric(alpha), beta, w, normalized, as.numeric(tol), max_iter)
}

hits_centrality <- function(graph, weights = NULL, tol = 1e-6, max_iter = 1000L) {
  pointer <- graph_pointer(graph)
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol <= 0) {
    stop("tol must be a positive number", call. = FALSE)
  }
  max_iter <- validate_count(max_iter)
  w <- resolve_edge_weights_opt(graph, weights)
  cpp_hits(pointer, w, as.numeric(tol), max_iter)
}

# ===========================================================================
# 2. Shortest Paths & Diameter
# ===========================================================================

shortest_distance <- function(graph, source = NULL, target = NULL,
                              weights = NULL, max_dist = Inf) {
  pointer <- graph_pointer(graph)
  nv <- as.integer(graph_num_vertices(graph))
  if (!is.numeric(max_dist) || length(max_dist) != 1L || is.na(max_dist) || max_dist < 0) {
    stop("max_dist must be a non-negative number (or Inf)", call. = FALSE)
  }
  srcs <- if (is.null(source)) {
    seq_len(nv)
  } else {
    if (!is.numeric(source) || anyNA(source) || any(source < 1 | source > nv | source != floor(source))) {
      stop("source must contain valid 1-based vertex IDs", call. = FALSE)
    }
    as.integer(source)
  }

  w <- resolve_edge_weights_opt(graph, weights)
  mat <- cpp_shortest_distance(pointer, srcs, w, as.numeric(max_dist))

  if (!is.null(target)) {
    if (!is.numeric(target) || anyNA(target) || any(target < 1 | target > nv | target != floor(target))) {
      stop("target must contain valid 1-based vertex IDs", call. = FALSE)
    }
    tgts <- as.integer(target)
    mat <- mat[, tgts, drop = FALSE]
  }
  if (nrow(mat) == 1L && !is.null(source) && length(source) == 1L) {
    return(as.numeric(mat[1L, ]))
  }
  mat
}

shortest_path <- function(graph, source, target, weights = NULL) {
  pointer <- graph_pointer(graph)
  nv <- as.integer(graph_num_vertices(graph))
  if (!is.numeric(source) || length(source) != 1L || is.na(source) ||
      source < 1 || source > nv || source != floor(source)) {
    stop("source must be a single valid 1-based vertex ID", call. = FALSE)
  }
  if (!is.numeric(target) || length(target) != 1L || is.na(target) ||
      target < 1 || target > nv || target != floor(target)) {
    stop("target must be a single valid 1-based vertex ID", call. = FALSE)
  }
  w <- resolve_edge_weights_opt(graph, weights)
  cpp_shortest_path(pointer, as.integer(source), as.integer(target), w)
}

pseudo_diameter <- function(graph, source = 1L, weights = NULL) {
  pointer <- graph_pointer(graph)
  nv <- as.integer(graph_num_vertices(graph))
  if (nv == 0L) return(list(diameter = 0, endpoints = integer(0)))
  if (!is.numeric(source) || length(source) != 1L || is.na(source) ||
      source < 1 || source > nv || source != floor(source)) {
    stop("source must be a single valid 1-based vertex ID", call. = FALSE)
  }
  w <- resolve_edge_weights_opt(graph, weights)
  cpp_pseudo_diameter(pointer, as.integer(source), w)
}

graph_diameter <- function(graph, weights = NULL) {
  dmat <- shortest_distance(graph, weights = weights)
  fin <- dmat[is.finite(dmat)]
  if (!length(fin)) return(0)
  max(fin)
}

# ===========================================================================
# 3. Connectivity, Components, DAG, Clustering & Assortativity
# ===========================================================================

label_components <- function(graph) {
  cpp_label_components(graph_pointer(graph))
}

label_biconnected_components <- function(graph) {
  cpp_biconnected_components(graph_pointer(graph))
}

topological_sort <- function(graph) {
  res <- cpp_topological_sort(graph_pointer(graph))
  if (!isTRUE(res$is_dag)) {
    stop("Graph contains a directed cycle and cannot be topologically sorted", call. = FALSE)
  }
  res$order
}

is_dag <- function(graph) {
  if (!isTRUE(graph$directed)) return(FALSE)
  isTRUE(cpp_topological_sort(graph_pointer(graph))$is_dag)
}

local_clustering <- function(graph) {
  cpp_clustering(graph_pointer(graph))$local
}

global_clustering <- function(graph) {
  res <- cpp_clustering(graph_pointer(graph))
  c(clustering = res$global, std_err = res$std_err)
}

assortativity <- function(graph, values = "total", categorical = FALSE, weights = NULL) {
  pointer <- graph_pointer(graph)
  nv <- as.integer(graph_num_vertices(graph))
  if (!is.logical(categorical) || length(categorical) != 1L || is.na(categorical)) {
    stop("categorical must be TRUE or FALSE", call. = FALSE)
  }

  if (is.character(values) && length(values) == 1L && values %in% c("total", "in", "out")) {
    if (!isTRUE(graph$directed)) {
      x_src <- x_tgt <- graph_degree(graph, "total")
    } else if (values == "total") {
      x_src <- x_tgt <- graph_degree(graph, "total")
    } else if (values == "out") {
      x_src <- x_tgt <- graph_degree(graph, "out")
    } else {
      x_src <- x_tgt <- graph_degree(graph, "in")
    }
  } else if (is.character(values) && length(values) == 1L &&
             values %in% graph_vertex_properties(graph)) {
    vp <- graph_vertex_property(graph, values)
    if (is.character(vp) || is.factor(vp) || is.logical(vp)) {
      categorical <- TRUE
      x_src <- x_tgt <- as.numeric(as.integer(as.factor(vp)))
    } else if (is.numeric(vp) && length(vp) == nv) {
      x_src <- x_tgt <- as.numeric(vp)
    } else {
      stop("vertex property for assortativity must be scalar atomic per vertex", call. = FALSE)
    }
  } else if (is.factor(values) || is.character(values) || is.logical(values)) {
    if (length(values) != nv) stop("values must have length graph_num_vertices(graph)", call. = FALSE)
    categorical <- TRUE
    x_src <- x_tgt <- as.numeric(as.integer(as.factor(values)))
  } else if (is.numeric(values) && length(values) == nv) {
    x_src <- x_tgt <- as.numeric(values)
  } else {
    stop("values must be 'total'/'in'/'out', a vertex property name, or a length-N vector",
         call. = FALSE)
  }

  w <- resolve_edge_weights_opt(graph, weights)
  cpp_assortativity(pointer, x_src, x_tgt, categorical, w)
}

# ===========================================================================
# 4. Classic Random Graph Generators
# ===========================================================================

random_graph_erdos_renyi <- function(n, p = NULL, m = NULL,
                                     directed = FALSE, self_loops = FALSE) {
  n <- validate_count(n)
  if (!is.logical(directed) || length(directed) != 1L || is.na(directed)) {
    stop("directed must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(self_loops) || length(self_loops) != 1L || is.na(self_loops)) {
    stop("self_loops must be TRUE or FALSE", call. = FALSE)
  }
  if (is.null(p) == is.null(m)) {
    stop("Specify exactly one of p (G(n, p)) or m (G(n, m))", call. = FALSE)
  }
  g <- graph_create(n, directed = directed)
  el <- if (!is.null(p)) {
    if (!is.numeric(p) || length(p) != 1L || is.na(p) || p < 0 || p > 1) {
      stop("p must be a probability in [0, 1]", call. = FALSE)
    }
    cpp_erdos_renyi_gnp(n, as.numeric(p), directed, self_loops)
  } else {
    m <- validate_count(m)
    cpp_erdos_renyi_gnm(n, m, directed, self_loops)
  }
  if (nrow(el) > 0L) graph_add_edges(g, el)
  g
}

random_graph_barabasi_albert <- function(n, m = 2L, gamma = 1.0, directed = FALSE) {
  n <- validate_count(n)
  m <- validate_count(m)
  if (!is.numeric(gamma) || length(gamma) != 1L || is.na(gamma) || !is.finite(gamma)) {
    stop("gamma must be a finite number", call. = FALSE)
  }
  if (!is.logical(directed) || length(directed) != 1L || is.na(directed)) {
    stop("directed must be TRUE or FALSE", call. = FALSE)
  }
  g <- graph_create(n, directed = directed)
  el <- cpp_barabasi_albert(n, m, as.numeric(gamma), directed)
  if (nrow(el) > 0L) graph_add_edges(g, el)
  g
}

random_graph_watts_strogatz <- function(n, k = 4L, p = 0.05, circular = TRUE) {
  n <- validate_count(n)
  k <- validate_count(k)
  if (!is.numeric(p) || length(p) != 1L || is.na(p) || p < 0 || p > 1) {
    stop("p must be a rewiring probability in [0, 1]", call. = FALSE)
  }
  if (!is.logical(circular) || length(circular) != 1L || is.na(circular)) {
    stop("circular must be TRUE or FALSE", call. = FALSE)
  }
  g <- graph_create(n, directed = FALSE)
  el <- cpp_watts_strogatz(n, k, as.numeric(p), circular)
  if (nrow(el) > 0L) graph_add_edges(g, el)
  g
}

random_graph_configuration <- function(out_deg, in_deg = NULL) {
  if (!is.numeric(out_deg) || anyNA(out_deg) || any(out_deg < 0 | out_deg != floor(out_deg))) {
    stop("out_deg must be a vector of non-negative integers", call. = FALSE)
  }
  n <- length(out_deg)
  directed <- !is.null(in_deg)
  if (!directed) {
    if (sum(out_deg) %% 2 != 0) {
      stop("Sum of undirected degrees in out_deg must be even", call. = FALSE)
    }
    g <- graph_create(n, directed = FALSE)
    stubs <- rep.int(seq_len(n), as.integer(out_deg))
    if (length(stubs) > 0L) {
      stubs <- sample(stubs)
      el <- matrix(stubs, ncol = 2L, byrow = TRUE)
      graph_add_edges(g, el)
    }
    return(g)
  }

  if (!is.numeric(in_deg) || length(in_deg) != n || anyNA(in_deg) ||
      any(in_deg < 0 | in_deg != floor(in_deg))) {
    stop("in_deg must be a non-negative integer vector of length(out_deg)", call. = FALSE)
  }
  if (sum(out_deg) != sum(in_deg)) {
    stop("sum(out_deg) must equal sum(in_deg) for directed configuration model", call. = FALSE)
  }
  g <- graph_create(n, directed = TRUE)
  s_out <- rep.int(seq_len(n), as.integer(out_deg))
  s_in <- sample(rep.int(seq_len(n), as.integer(in_deg)))
  if (length(s_out) > 0L) {
    graph_add_edges(g, cbind(s_out, s_in))
  }
  g
}

random_graph_sbm <- function(blocks, probs, directed = FALSE, self_loops = FALSE) {
  if (!is.numeric(blocks) || !length(blocks) || anyNA(blocks) ||
      any(blocks < 1 | blocks != floor(blocks))) {
    stop("blocks must be a non-empty vector of positive integer group labels", call. = FALSE)
  }
  b <- match(as.integer(blocks), sort(unique(as.integer(blocks))))
  B <- max(b)
  if (!is.matrix(probs) || !identical(dim(probs), c(B, B)) ||
      anyNA(probs) || any(probs < 0 | probs > 1)) {
    stop(sprintf("probs must be a %d x %d probability matrix in [0, 1]", B, B),
         call. = FALSE)
  }
  if (!directed && !isSymmetric(unname(probs), tol = 1e-12)) {
    stop("probs must be symmetric when directed = FALSE", call. = FALSE)
  }

  n <- length(b)
  g <- graph_create(n, directed = directed)
  by_blk <- split(seq_len(n), b)
  el_parts <- list()
  idx <- 1L
  for (r in seq_len(B)) {
    vr <- by_blk[[r]]
    s_range <- if (directed) seq_len(B) else r:B
    for (s in s_range) {
      p <- probs[r, s]
      if (p <= 0) next
      vs <- by_blk[[s]]
      if (r == s) {
        sub <- cpp_erdos_renyi_gnp(length(vr), p, directed, self_loops)
        if (nrow(sub) > 0L) {
          el_parts[[idx]] <- cbind(vr[sub[, 1L]], vr[sub[, 2L]])
          idx <- idx + 1L
        }
      } else {
        # Rectangular Bernoulli block vr x vs
        m_expected <- length(vr) * length(vs)
        draw <- which(matrix(stats::runif(m_expected) < p,
                             nrow = length(vr), ncol = length(vs)),
                      arr.ind = TRUE)
        if (nrow(draw) > 0L) {
          el_parts[[idx]] <- cbind(vr[draw[, 1L]], vs[draw[, 2L]])
          idx <- idx + 1L
        }
      }
    }
  }
  if (length(el_parts)) {
    graph_add_edges(g, do.call(rbind, el_parts))
  }
  graph_set_vertex_property(g, "block", b)
  g
}

# ===========================================================================
# 5. Modularity & Modularity Maximization (WITH EXPLICIT STATISTICAL WARNINGS)
# ===========================================================================

emit_modularity_warning <- function() {
  warning(
    "WARNING -- Modularity maximization is almost always a terrible idea!\n",
    "  Modularity maximization lacks any statistical regularization: it massively\n",
    "  overfits sparse and random networks (finding spurious high-modularity\n",
    "  'communities' in pure Erdos-Renyi and configuration-model graphs) while\n",
    "  simultaneously underfitting real small communities due to the Fortunato-\n",
    "  Barthelemy resolution limit.\n",
    "  It is included in rgraphtool strictly as a baseline for comparison.\n",
    "  For principled community detection, use minimize_blockmodel_dl(),\n",
    "  minimize_nested_blockmodel_dl(), or minimize_pp_blockmodel_dl() instead.\n",
    "  See Tiago P. Peixoto (2023), 'Descriptive vs. inferential community\n",
    "  detection in networks: pitfalls, myths and half-truths' and\n",
    "  https://skewed.de/tiago/blog/modularity-harmful (suppress with warn = FALSE).",
    call. = FALSE
  )
}

modularity <- function(graph, blocks, gamma = 1.0, weights = NULL) {
  pointer <- graph_pointer(graph)
  nv <- as.integer(graph_num_vertices(graph))
  if (inherits(blocks, "rgraphtool_modularity_state")) {
    blocks <- blocks$blocks
  } else if (inherits(blocks, "rgraphtool_block_state") ||
             inherits(blocks, "rgraphtool_pp_block_state") ||
             inherits(blocks, "rgraphtool_weighted_block_state")) {
    blocks <- block_membership(blocks)
  }
  if (!is.numeric(blocks) || length(blocks) != nv || anyNA(blocks) ||
      any(!is.finite(blocks)) || any(blocks != floor(blocks))) {
    stop("blocks must be an integer partition vector of length graph_num_vertices(graph)",
         call. = FALSE)
  }
  if (!is.numeric(gamma) || length(gamma) != 1L || is.na(gamma) || !is.finite(gamma)) {
    stop("gamma must be a finite scalar resolution parameter", call. = FALSE)
  }
  w <- resolve_edge_weights_opt(graph, weights)
  cpp_modularity(pointer, as.integer(blocks), as.numeric(gamma), w)
}

maximize_modularity <- function(graph, gamma = 1.0, weights = NULL,
                                nruns = 1L, max_levels = 10L, warn = TRUE) {
  pointer <- graph_pointer(graph)
  if (!is.numeric(gamma) || length(gamma) != 1L || is.na(gamma) || !is.finite(gamma)) {
    stop("gamma must be a finite scalar resolution parameter", call. = FALSE)
  }
  nruns <- as.integer(nruns)
  max_levels <- as.integer(max_levels)
  if (is.na(nruns) || nruns < 1L) stop("nruns must be a positive integer", call. = FALSE)
  if (is.na(max_levels) || max_levels < 1L) stop("max_levels must be a positive integer", call. = FALSE)
  if (isTRUE(warn)) {
    emit_modularity_warning()
  }
  w <- resolve_edge_weights_opt(graph, weights)
  res <- cpp_maximize_modularity(pointer, as.numeric(gamma), w, nruns, max_levels)
  structure(list(
    blocks = res$blocks,
    modularity = res$modularity,
    B = res$B,
    gamma = as.numeric(gamma),
    graph = graph
  ), class = "rgraphtool_modularity_state")
}

print.rgraphtool_modularity_state <- function(x, ...) {
  cat(sprintf("<modularity maximization state (gamma = %g): %d groups, Q = %.4f>\n",
              x$gamma, x$B, x$modularity))
  cat("  WARNING: Modularity maximization lacks statistical regularization and overfits\n")
  cat("  random graphs. Prefer minimize_blockmodel_dl() or minimize_pp_blockmodel_dl().\n")
  invisible(x)
}
