# Weight-informed starting partitions for the weighted minimizers.
#
# The agglomerative minimizers propose merges from neighbours' groups, which
# carry no information when only the edge weights are structured (e.g. a dense
# or random topology with block-dependent weights): the search then stops at
# B = 1 although structured partitions have a lower description length. These
# seeds give the minimizer extra starts that see the weights; it refines each
# and keeps the lowest description length, so they can only help.
#
# Embedding: leading eigenvectors (largest |eigenvalue|) of the weighted
# adjacency with the weight centred over edges (log1p for non-negative types),
# by subspace iteration over the edge list (no N x N matrix). Seeds for each B:
# varimax rotation of the top-B eigenvectors (Rohe & Zeng 2023: columns made
# positively skewed, each vertex to its largest loading) and k-means.

.weight_transform <- function(w, type) {
  if (type == "real-normal") return(w)
  if (any(w < 0)) return(w)
  log1p(w)
}

.edge_matvec <- function(X, u, v, x, n) {
  Y <- matrix(0, n, ncol(X))
  a <- rowsum(x * X[v, , drop = FALSE], u, reorder = FALSE)
  Y[as.integer(rownames(a)), ] <- Y[as.integer(rownames(a)), , drop = FALSE] + a
  off <- u != v
  if (any(off)) {
    b <- rowsum(x[off] * X[u[off], , drop = FALSE], v[off], reorder = FALSE)
    Y[as.integer(rownames(b)), ] <- Y[as.integer(rownames(b)), , drop = FALSE] + b
  }
  Y
}

.spectral_embedding <- function(u, v, x, n, p, iters = 40L) {
  X <- qr.Q(qr(matrix(stats::rnorm(n * p), n, p)))
  for (it in seq_len(iters)) X <- qr.Q(qr(.edge_matvec(X, u, v, x, n)))
  AX <- .edge_matvec(X, u, v, x, n)
  ev <- eigen(crossprod(X, AX), symmetric = TRUE)
  U <- X %*% ev$vectors
  U[, order(abs(ev$values), decreasing = TRUE), drop = FALSE]
}

.varimax_seed <- function(U, B) {
  Ub <- sqrt(nrow(U)) * U[, seq_len(B), drop = FALSE]
  R <- unclass(stats::varimax(Ub, normalize = FALSE)$loadings)
  R <- sweep(R, 2L, ifelse(colSums(R^3) < 0, -1, 1), "*")
  max.col(R, ties.method = "first")
}

weight_spectral_seeds <- function(graph, wmat, types, B_max = 8L) {
  # Leave the caller's random stream untouched: the minimizer's own starts must
  # be the same with and without these seeds, so the extra starts can only
  # lower the description length.
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
  on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
          else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
            rm(".Random.seed", envir = globalenv()), add = TRUE)
  n <- graph_num_vertices(graph)
  if (n < 10L || nrow(wmat) < n || ncol(wmat) < 1L) return(NULL)
  ed <- graph_edges(graph)
  x <- .weight_transform(wmat[, 1L], types[1L])
  x <- x - mean(x)
  if (!any(abs(x) > 0)) return(NULL)
  B_max <- as.integer(min(B_max, n - 1L))
  if (B_max < 2L) return(NULL)
  U <- .spectral_embedding(as.integer(ed[, 1L]), as.integer(ed[, 2L]), x, n,
                           p = min(n - 1L, B_max + 4L))
  seeds <- list()
  for (B in 2:B_max) {
    seeds[[length(seeds) + 1L]] <- .varimax_seed(U, B)
    km <- tryCatch(suppressWarnings(stats::kmeans(U[, seq_len(B), drop = FALSE], B, nstart = 5L,
                                                   iter.max = 50L))$cluster,
                   error = function(e) NULL)
    if (!is.null(km)) seeds[[length(seeds) + 1L]] <- km
  }
  seeds <- lapply(seeds, function(b) as.integer(match(b, unique(b))))
  seeds[!duplicated(vapply(seeds, paste, "", collapse = ","))]
}
