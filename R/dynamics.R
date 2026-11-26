# Network reconstruction from dynamics and behavior (Peixoto 2019 PRL,
# Peixoto 2024, Peixoto 2025 PRX), L1 regularization, and forward simulators.
#
# Copyright (C) 2006-2026 Tiago de Paula Peixoto (graph-tool algorithms)
# Copyright (C) 2026 Diogo Melo (R port)
# LGPL-3+

.DYN_PSEUDO_ISING   <- 0L
.DYN_ISING_GLAUBER  <- 1L
.DYN_PSEUDO_CISING  <- 2L
.DYN_CISING_GLAUBER <- 3L
.DYN_EPIDEMICS      <- 4L
.DYN_PSEUDO_NORMAL  <- 5L
.DYN_NORMAL_GLAUBER <- 6L
.DYN_LINEAR_NORMAL  <- 7L
.DYN_PSEUDO_POTTS   <- 8L
.DYN_POTTS_GLAUBER  <- 9L

.dyn_normalize_s <- function(s, fix_mean = FALSE) {
  if (!is.null(attr(s, "gram"))) {
    return(s)
  }
  if (is.matrix(s)) {
    s_list <- list(s)
  } else if (is.list(s) && length(s) >= 1L) {
    s_list <- s
  } else {
    stop("s must be an N x M numeric/integer matrix or a non-empty list of N x M_c matrices",
         call. = FALSE)
  }
  N <- NULL
  out <- vector("list", length(s_list))
  for (idx in seq_along(s_list)) {
    m <- s_list[[idx]]
    if (is.data.frame(m)) m <- as.matrix(m)
    if (!is.matrix(m) || !is.numeric(m) || anyNA(m) || any(!is.finite(m))) {
      stop("Each element of s must be a finite numeric matrix without NAs", call. = FALSE)
    }
    if (is.null(N)) N <- nrow(m)
    if (nrow(m) != N) stop("All matrices in s must have the same number of rows N", call. = FALSE)
    if (ncol(m) < 1L) stop("Data matrices in s must have at least one column", call. = FALSE)
    storage.mode(m) <- "double"
    out[[idx]] <- m
  }
  if (isTRUE(fix_mean)) {
    # Center each vertex across all observations
    tot_sum <- numeric(N)
    tot_cnt <- 0L
    for (m in out) {
      tot_sum <- tot_sum + rowSums(m)
      tot_cnt <- tot_cnt + ncol(m)
    }
    mu <- tot_sum / tot_cnt
    out <- lapply(out, function(m) m - mu)
  }
  out
}

.dyn_resolve_theta <- function(theta, N, q_t) {
  if (is.null(theta)) {
    return(matrix(0, nrow = N, ncol = q_t))
  }
  if (is.matrix(theta)) {
    if (nrow(theta) != N || ncol(theta) != q_t || !is.numeric(theta) || anyNA(theta)) {
      stop(sprintf("theta matrix must be numeric of dimension %d x %d", N, q_t), call. = FALSE)
    }
    storage.mode(theta) <- "double"
    return(theta)
  }
  if (is.numeric(theta) && !anyNA(theta)) {
    if (length(theta) == 1L) {
      return(matrix(as.numeric(theta), nrow = N, ncol = q_t))
    }
    if (length(theta) == q_t) {
      return(matrix(rep(as.numeric(theta), each = N), nrow = N, ncol = q_t))
    }
    if (q_t == 1L && length(theta) == N) {
      return(matrix(as.numeric(theta), nrow = N, ncol = 1L))
    }
  }
  stop(sprintf("theta must be NULL, length 1, length %d, or an %d x %d numeric matrix",
               q_t, N, q_t), call. = FALSE)
}

.dyn_resolve_range_list <- function(rng, q_t, name) {
  if (is.numeric(rng) && length(rng) == 2L) {
    return(list(lo = rep(as.numeric(rng[1L]), q_t), hi = rep(as.numeric(rng[2L]), q_t)))
  }
  if (is.list(rng) && length(rng) == q_t) {
    lo <- vapply(rng, function(z) as.numeric(z[1L]), numeric(1L))
    hi <- vapply(rng, function(z) as.numeric(z[2L]), numeric(1L))
    return(list(lo = lo, hi = hi))
  }
  stop(sprintf("%s must be a length-2 numeric vector or a list of %d length-2 bounds", name, q_t),
       call. = FALSE)
}

.dyn_build_state <- function(model_code, model_name, s, graph, directed, self_loops,
                             q_t, f_mat, has_zero, exposed, positive, pslack,
                             x, x_range, theta, theta_range,
                             xdelta, tdelta, disable_xdist, disable_tdist,
                             xl1, tl1, sbm, nested, blocks,
                             degree_corrected, degree_prior, fix_mean, q_max) {
  s_list <- .dyn_normalize_s(s, fix_mean = fix_mean)
  gram_attr <- attr(s_list, "gram")
  N <- nrow(s_list[[1L]])
  n_obs_val <- if (!is.null(gram_attr)) {
    as.numeric(gram_attr$n_eff)
  } else {
    sum(vapply(s_list, ncol, integer(1L)))
  }
  degree_prior <- match.arg(degree_prior, c("distributed", "uniform"))
  dist_flag <- identical(degree_prior, "distributed")
  nested <- isTRUE(nested)
  degree_corrected <- isTRUE(degree_corrected)
  self_loops <- isTRUE(self_loops)

  if (!is.null(graph)) {
    graph_pointer(graph)
    if (graph_num_vertices(graph) != N) {
      stop(sprintf("graph has %d vertices, expected N = %d from s", graph_num_vertices(graph), N),
           call. = FALSE)
    }
    directed <- isTRUE(graph$directed)
    init_e <- graph_edges(graph)
    if (nrow(init_e) == 0L) {
      init_x <- numeric(0)
    } else if (is.null(x)) {
      init_x <- rep(if (model_code == .DYN_EPIDEMICS) -0.5 else 1.0, nrow(init_e))
    } else if (is.character(x) && length(x) == 1L) {
      init_x <- as.numeric(graph_edge_property(graph, x))
    } else if (is.numeric(x) && length(x) == 1L) {
      init_x <- rep(as.numeric(x), nrow(init_e))
    } else if (is.numeric(x) && length(x) == nrow(init_e)) {
      init_x <- as.numeric(x)
    } else {
      stop("x must be NULL, scalar, edge property name, or numeric vector matching graph edges",
           call. = FALSE)
    }
  } else {
    directed <- isTRUE(directed)
    init_e <- matrix(integer(0), nrow = 0L, ncol = 2L)
    init_x <- numeric(0)
  }

  th_mat <- .dyn_resolve_theta(theta, N, q_t)
  tr <- .dyn_resolve_range_list(theta_range, q_t, "theta_range")
  if (!is.numeric(x_range) || length(x_range) != 2L || x_range[1L] >= x_range[2L]) {
    stop("x_range must be a length-2 numeric vector with x_range[1] < x_range[2]", call. = FALSE)
  }
  bs <- .uncertain_validate_bs(N, blocks, nested)

  tdelta_v <- if (length(tdelta) == 1L) rep(as.numeric(tdelta), q_t) else as.numeric(tdelta)
  tl1_v <- if (length(tl1) == 1L) rep(as.numeric(tl1), q_t) else as.numeric(tl1)
  if (length(tdelta_v) != q_t || length(tl1_v) != q_t) {
    stop(sprintf("tdelta and tl1 must have length 1 or %d", q_t), call. = FALSE)
  }

  xdist_def <- !isTRUE(disable_xdist)
  tdist_def <- rep(!isTRUE(disable_tdist), q_t)
  sbm_def <- if (is.null(sbm)) !isTRUE(disable_xdist) else isTRUE(sbm)

  structure(
    list(
      model_code = as.integer(model_code),
      model = model_name,
      n_vertices = as.integer(N),
      directed = directed,
      self_loops = self_loops,
      q_t = as.integer(q_t),
      s_list = s_list,
      n_obs = n_obs_val,
      f_mat = f_mat,
      has_zero = isTRUE(has_zero),
      exposed = isTRUE(exposed),
      positive = isTRUE(positive),
      pslack = as.numeric(pslack),
      x_range = as.numeric(x_range),
      tmin = tr$lo,
      tmax = tr$hi,
      edges = init_e,
      x = init_x,
      theta = th_mat,
      bs = bs,
      nested = nested,
      degree_corrected = degree_corrected,
      degree_prior = degree_prior,
      dist = dist_flag,
      q_max = as.numeric(q_max),
      default_prior = list(
        sbm = sbm_def,
        latent_edges = TRUE,
        density = FALSE,
        aE = 1.0,
        xdist = xdist_def,
        xdist_uniform = FALSE,
        xl1 = as.numeric(xl1),
        xdelta = as.numeric(xdelta),
        tdist = tdist_def,
        tdist_uniform = rep(FALSE, q_t),
        tl1 = tl1_v,
        tdelta = tdelta_v,
        normal = FALSE,
        mu = 0.0,
        sigma = 1.0
      )
    ),
    class = "rgraphtool_dynamics_block_state"
  )
}

# ---- Constructors ----

#' @noRd
#' @export
pseudo_ising_block_state <- function(s, graph = NULL, has_zero = FALSE,
                                     x = 1, x_range = c(-Inf, Inf),
                                     theta = 0, theta_range = c(-Inf, Inf),
                                     xdelta = 1e-8, tdelta = 1e-8,
                                     disable_xdist = FALSE, disable_tdist = FALSE,
                                     xl1 = 1, tl1 = 1, sbm = NULL,
                                     nested = TRUE, blocks = NULL,
                                     degree_corrected = TRUE,
                                     degree_prior = c("distributed", "uniform"),
                                     self_loops = FALSE, q_max = Inf) {
  .dyn_build_state(
    .DYN_PSEUDO_ISING, "pseudo_ising", s, graph, directed = FALSE, self_loops = self_loops,
    q_t = 1L, f_mat = matrix(0, 0, 0), has_zero = has_zero, exposed = FALSE,
    positive = FALSE, pslack = 1e-6, x = x, x_range = x_range,
    theta = theta, theta_range = theta_range, xdelta = xdelta, tdelta = tdelta,
    disable_xdist = disable_xdist, disable_tdist = disable_tdist,
    xl1 = xl1, tl1 = tl1, sbm = sbm, nested = nested, blocks = blocks,
    degree_corrected = degree_corrected, degree_prior = match.arg(degree_prior),
    fix_mean = FALSE, q_max = q_max
  )
}

#' @noRd
#' @export
ising_glauber_block_state <- function(s, graph = NULL, directed = TRUE, has_zero = FALSE,
                                      x = 1, x_range = c(-Inf, Inf),
                                      theta = 0, theta_range = c(-Inf, Inf),
                                      xdelta = 1e-8, tdelta = 1e-8,
                                      disable_xdist = FALSE, disable_tdist = FALSE,
                                      xl1 = 1, tl1 = 1, sbm = NULL,
                                      nested = TRUE, blocks = NULL,
                                      degree_corrected = TRUE,
                                      degree_prior = c("distributed", "uniform"),
                                      self_loops = FALSE, q_max = Inf) {
  .dyn_build_state(
    .DYN_ISING_GLAUBER, "ising_glauber", s, graph, directed = directed, self_loops = self_loops,
    q_t = 1L, f_mat = matrix(0, 0, 0), has_zero = has_zero, exposed = FALSE,
    positive = FALSE, pslack = 1e-6, x = x, x_range = x_range,
    theta = theta, theta_range = theta_range, xdelta = xdelta, tdelta = tdelta,
    disable_xdist = disable_xdist, disable_tdist = disable_tdist,
    xl1 = xl1, tl1 = tl1, sbm = sbm, nested = nested, blocks = blocks,
    degree_corrected = degree_corrected, degree_prior = match.arg(degree_prior),
    fix_mean = FALSE, q_max = q_max
  )
}

#' @noRd
#' @export
pseudo_cising_block_state <- function(s, graph = NULL,
                                      x = 1, x_range = c(-Inf, Inf),
                                      theta = 0, theta_range = c(-Inf, Inf),
                                      xdelta = 1e-8, tdelta = 1e-8,
                                      disable_xdist = FALSE, disable_tdist = FALSE,
                                      xl1 = 1, tl1 = 1, sbm = NULL,
                                      nested = TRUE, blocks = NULL,
                                      degree_corrected = TRUE,
                                      degree_prior = c("distributed", "uniform"),
                                      self_loops = FALSE, q_max = Inf) {
  .dyn_build_state(
    .DYN_PSEUDO_CISING, "pseudo_cising", s, graph, directed = FALSE, self_loops = self_loops,
    q_t = 1L, f_mat = matrix(0, 0, 0), has_zero = FALSE, exposed = FALSE,
    positive = FALSE, pslack = 1e-6, x = x, x_range = x_range,
    theta = theta, theta_range = theta_range, xdelta = xdelta, tdelta = tdelta,
    disable_xdist = disable_xdist, disable_tdist = disable_tdist,
    xl1 = xl1, tl1 = tl1, sbm = sbm, nested = nested, blocks = blocks,
    degree_corrected = degree_corrected, degree_prior = match.arg(degree_prior),
    fix_mean = FALSE, q_max = q_max
  )
}

#' @noRd
#' @export
cising_glauber_block_state <- function(s, graph = NULL, directed = TRUE,
                                       x = 1, x_range = c(-Inf, Inf),
                                       theta = 0, theta_range = c(-Inf, Inf),
                                       xdelta = 1e-8, tdelta = 1e-8,
                                       disable_xdist = FALSE, disable_tdist = FALSE,
                                       xl1 = 1, tl1 = 1, sbm = NULL,
                                       nested = TRUE, blocks = NULL,
                                       degree_corrected = TRUE,
                                       degree_prior = c("distributed", "uniform"),
                                       self_loops = FALSE, q_max = Inf) {
  .dyn_build_state(
    .DYN_CISING_GLAUBER, "cising_glauber", s, graph, directed = directed, self_loops = self_loops,
    q_t = 1L, f_mat = matrix(0, 0, 0), has_zero = FALSE, exposed = FALSE,
    positive = FALSE, pslack = 1e-6, x = x, x_range = x_range,
    theta = theta, theta_range = theta_range, xdelta = xdelta, tdelta = tdelta,
    disable_xdist = disable_xdist, disable_tdist = disable_tdist,
    xl1 = xl1, tl1 = tl1, sbm = sbm, nested = nested, blocks = blocks,
    degree_corrected = degree_corrected, degree_prior = match.arg(degree_prior),
    fix_mean = FALSE, q_max = q_max
  )
}

#' @noRd
#' @export
epidemics_block_state <- function(s, graph = NULL, directed = TRUE, exposed = FALSE,
                                  x = -0.5, x_range = c(-Inf, -1e-8),
                                  theta = -3, theta_range = c(-Inf, -1e-8),
                                  xdelta = 1e-8, tdelta = 1e-8,
                                  disable_xdist = FALSE, disable_tdist = FALSE,
                                  xl1 = 1, tl1 = 1, sbm = NULL,
                                  nested = TRUE, blocks = NULL,
                                  degree_corrected = TRUE,
                                  degree_prior = c("distributed", "uniform"),
                                  self_loops = FALSE, q_max = Inf) {
  .dyn_build_state(
    .DYN_EPIDEMICS, "epidemics", s, graph, directed = directed, self_loops = self_loops,
    q_t = 1L, f_mat = matrix(0, 0, 0), has_zero = FALSE, exposed = exposed,
    positive = FALSE, pslack = 1e-6, x = x, x_range = x_range,
    theta = theta, theta_range = theta_range, xdelta = xdelta, tdelta = tdelta,
    disable_xdist = disable_xdist, disable_tdist = disable_tdist,
    xl1 = xl1, tl1 = tl1, sbm = sbm, nested = nested, blocks = blocks,
    degree_corrected = degree_corrected, degree_prior = match.arg(degree_prior),
    fix_mean = FALSE, q_max = q_max
  )
}

#' @noRd
#' @export
pseudo_normal_block_state <- function(s = NULL, graph = NULL, fix_mean = FALSE,
                                      positive = TRUE, pslack = 1e-6,
                                      x = 0.1, x_range = c(-Inf, Inf),
                                      theta = c(0, 0),
                                      theta_range = list(c(-200, 200), c(-Inf, Inf)),
                                      xdelta = 1e-8, tdelta = c(1e-8, 1e-8),
                                      disable_xdist = FALSE, disable_tdist = FALSE,
                                      xl1 = 1, tl1 = c(1, 1), sbm = NULL,
                                      nested = TRUE, blocks = NULL,
                                      degree_corrected = TRUE,
                                      degree_prior = c("distributed", "uniform"),
                                      self_loops = FALSE, q_max = Inf,
                                      cov = NULL, n_eff = NULL, mean = NULL) {
  has_s <- !is.null(s)
  has_cov <- !is.null(cov)
  if (has_s == has_cov) {
    stop("Specify exactly one of 's' (sample matrix) or 'cov' (covariance matrix with 'n_eff')",
         call. = FALSE)
  }
  if (has_cov) {
    if (is.data.frame(cov)) cov <- as.matrix(cov)
    if (!is.matrix(cov) || !is.numeric(cov) || nrow(cov) != ncol(cov) ||
        nrow(cov) < 1L || anyNA(cov) || any(!is.finite(cov))) {
      stop("cov must be a finite square N x N numeric matrix", call. = FALSE)
    }
    if (!isSymmetric(unname(cov), tol = 1e-10)) {
      stop("cov must be symmetric", call. = FALSE)
    }
    if (any(diag(cov) <= 0)) {
      stop("cov must have strictly positive diagonal entries", call. = FALSE)
    }
    if (is.null(n_eff) || !is.numeric(n_eff) || length(n_eff) != 1L ||
        is.na(n_eff) || !is.finite(n_eff) || n_eff <= 0) {
      stop("n_eff must be a single positive finite number when cov is provided",
           call. = FALSE)
    }
    N_cov <- nrow(cov)
    if (is.null(mean)) {
      mean_vec <- numeric(N_cov)
    } else {
      if (!is.numeric(mean) || length(mean) != N_cov || anyNA(mean) || any(!is.finite(mean))) {
        stop("mean must be a finite numeric vector of length nrow(cov)", call. = FALSE)
      }
      mean_vec <- as.numeric(mean)
    }
    storage.mode(cov) <- "double"
    cov <- 0.5 * (cov + t(cov))
    s <- list(matrix(0, nrow = N_cov, ncol = 0L))
    attr(s, "gram") <- list(
      cov = cov,
      n_eff = as.numeric(n_eff),
      mean = if (isTRUE(fix_mean)) numeric(N_cov) else mean_vec
    )
  }
  .dyn_build_state(
    .DYN_PSEUDO_NORMAL, "pseudo_normal", s, graph, directed = FALSE, self_loops = self_loops,
    q_t = 2L, f_mat = matrix(0, 0, 0), has_zero = FALSE, exposed = FALSE,
    positive = positive, pslack = pslack, x = x, x_range = x_range,
    theta = theta, theta_range = theta_range, xdelta = xdelta, tdelta = tdelta,
    disable_xdist = disable_xdist, disable_tdist = disable_tdist,
    xl1 = xl1, tl1 = tl1, sbm = sbm, nested = nested, blocks = blocks,
    degree_corrected = degree_corrected, degree_prior = match.arg(degree_prior),
    fix_mean = fix_mean, q_max = q_max
  )
}

#' @noRd
#' @export
normal_glauber_block_state <- function(s, graph = NULL, directed = TRUE, fix_mean = FALSE,
                                       x = 0.1, x_range = c(-Inf, Inf),
                                       theta = c(0, 0),
                                       theta_range = list(c(-200, 200), c(-Inf, Inf)),
                                       xdelta = 1e-8, tdelta = c(1e-8, 1e-8),
                                       disable_xdist = FALSE, disable_tdist = FALSE,
                                       xl1 = 1, tl1 = c(1, 1), sbm = NULL,
                                       nested = TRUE, blocks = NULL,
                                       degree_corrected = TRUE,
                                       degree_prior = c("distributed", "uniform"),
                                       self_loops = FALSE, q_max = Inf) {
  .dyn_build_state(
    .DYN_NORMAL_GLAUBER, "normal_glauber", s, graph, directed = directed, self_loops = self_loops,
    q_t = 2L, f_mat = matrix(0, 0, 0), has_zero = FALSE, exposed = FALSE,
    positive = FALSE, pslack = 1e-6, x = x, x_range = x_range,
    theta = theta, theta_range = theta_range, xdelta = xdelta, tdelta = tdelta,
    disable_xdist = disable_xdist, disable_tdist = disable_tdist,
    xl1 = xl1, tl1 = tl1, sbm = sbm, nested = nested, blocks = blocks,
    degree_corrected = degree_corrected, degree_prior = match.arg(degree_prior),
    fix_mean = fix_mean, q_max = q_max
  )
}

#' @noRd
#' @export
linear_normal_block_state <- function(s, graph = NULL, directed = TRUE, self_loops = TRUE,
                                      x = 0.1, x_range = c(-Inf, Inf),
                                      theta = 0, theta_range = c(-200, 200),
                                      xdelta = 1e-8, tdelta = 1e-8,
                                      disable_xdist = FALSE, disable_tdist = FALSE,
                                      xl1 = 1, tl1 = 1, sbm = NULL,
                                      nested = TRUE, blocks = NULL,
                                      degree_corrected = TRUE,
                                      degree_prior = c("distributed", "uniform"),
                                      q_max = Inf) {
  .dyn_build_state(
    .DYN_LINEAR_NORMAL, "linear_normal", s, graph, directed = directed, self_loops = self_loops,
    q_t = 1L, f_mat = matrix(0, 0, 0), has_zero = FALSE, exposed = FALSE,
    positive = FALSE, pslack = 1e-6, x = x, x_range = x_range,
    theta = theta, theta_range = theta_range, xdelta = xdelta, tdelta = tdelta,
    disable_xdist = disable_xdist, disable_tdist = disable_tdist,
    xl1 = xl1, tl1 = tl1, sbm = sbm, nested = nested, blocks = blocks,
    degree_corrected = degree_corrected, degree_prior = match.arg(degree_prior),
    fix_mean = FALSE, q_max = q_max
  )
}

.dyn_resolve_potts_f <- function(s_list, f, q) {
  if (is.null(f)) {
    if (is.null(q)) {
      q <- max(vapply(s_list, function(m) as.integer(max(m)), integer(1L))) + 1L
    }
    q <- as.integer(q)
    if (is.na(q) || q < 2L) stop("Potts q must be an integer >= 2", call. = FALSE)
    f_mat <- diag(1, q)
  } else {
    f_mat <- as.matrix(f)
    if (!is.numeric(f_mat) || nrow(f_mat) != ncol(f_mat) || nrow(f_mat) < 2L || anyNA(f_mat)) {
      stop("f must be a square numeric matrix of dimension q x q (q >= 2)", call. = FALSE)
    }
    storage.mode(f_mat) <- "double"
  }
  q_out <- nrow(f_mat)
  for (m in s_list) {
    if (any(m < 0 | m >= q_out | m != floor(m))) {
      stop(sprintf("Potts states in s must be integers in 0..%d", q_out - 1L), call. = FALSE)
    }
  }
  f_mat
}

#' @noRd
#' @export
pseudo_potts_block_state <- function(s, f = NULL, q = NULL, graph = NULL,
                                     x = 1, x_range = c(-Inf, Inf),
                                     theta = NULL, theta_range = c(-Inf, Inf),
                                     xdelta = 1e-8, tdelta = 1e-8,
                                     disable_xdist = FALSE, disable_tdist = FALSE,
                                     xl1 = 1, tl1 = 1, sbm = NULL,
                                     nested = TRUE, blocks = NULL,
                                     degree_corrected = TRUE,
                                     degree_prior = c("distributed", "uniform"),
                                     self_loops = FALSE, q_max = Inf) {
  s_list <- .dyn_normalize_s(s, fix_mean = FALSE)
  f_mat <- .dyn_resolve_potts_f(s_list, f, q)
  q_t <- nrow(f_mat)
  .dyn_build_state(
    .DYN_PSEUDO_POTTS, "pseudo_potts", s_list, graph, directed = FALSE, self_loops = self_loops,
    q_t = q_t, f_mat = f_mat, has_zero = FALSE, exposed = FALSE,
    positive = FALSE, pslack = 1e-6, x = x, x_range = x_range,
    theta = theta, theta_range = theta_range, xdelta = xdelta, tdelta = tdelta,
    disable_xdist = disable_xdist, disable_tdist = disable_tdist,
    xl1 = xl1, tl1 = tl1, sbm = sbm, nested = nested, blocks = blocks,
    degree_corrected = degree_corrected, degree_prior = match.arg(degree_prior),
    fix_mean = FALSE, q_max = q_max
  )
}

#' @noRd
#' @export
potts_glauber_block_state <- function(s, f = NULL, q = NULL, graph = NULL, directed = TRUE,
                                      x = 1, x_range = c(-Inf, Inf),
                                      theta = NULL, theta_range = c(-Inf, Inf),
                                      xdelta = 1e-8, tdelta = 1e-8,
                                      disable_xdist = FALSE, disable_tdist = FALSE,
                                      xl1 = 1, tl1 = 1, sbm = NULL,
                                      nested = TRUE, blocks = NULL,
                                      degree_corrected = TRUE,
                                      degree_prior = c("distributed", "uniform"),
                                      self_loops = FALSE, q_max = Inf) {
  s_list <- .dyn_normalize_s(s, fix_mean = FALSE)
  f_mat <- .dyn_resolve_potts_f(s_list, f, q)
  q_t <- nrow(f_mat)
  .dyn_build_state(
    .DYN_POTTS_GLAUBER, "potts_glauber", s_list, graph, directed = directed, self_loops = self_loops,
    q_t = q_t, f_mat = f_mat, has_zero = FALSE, exposed = FALSE,
    positive = FALSE, pslack = 1e-6, x = x, x_range = x_range,
    theta = theta, theta_range = theta_range, xdelta = xdelta, tdelta = tdelta,
    disable_xdist = disable_xdist, disable_tdist = disable_tdist,
    xl1 = xl1, tl1 = tl1, sbm = sbm, nested = nested, blocks = blocks,
    degree_corrected = degree_corrected, degree_prior = match.arg(degree_prior),
    fix_mean = FALSE, q_max = q_max
  )
}

# ---- Evaluation, reconstruction, and L1 regularization path ----

.check_dynamics_state <- function(state) {
  if (!inherits(state, "rgraphtool_dynamics_block_state")) {
    stop("state must be created by a dynamics block state constructor (e.g., pseudo_ising_block_state())",
         call. = FALSE)
  }
}

.dyn_merge_prior <- function(state, sbm = NULL, latent_edges = NULL, density = NULL, aE = NULL,
                             xdist = NULL, tdist = NULL, xdist_uniform = NULL,
                             tdist_uniform = NULL, xl1 = NULL, tl1 = NULL,
                             normal = NULL, mu = NULL, sigma = NULL) {
  p <- state$default_prior
  if (!is.null(sbm)) p$sbm <- isTRUE(sbm)
  if (!is.null(latent_edges)) p$latent_edges <- isTRUE(latent_edges)
  if (!is.null(density)) p$density <- isTRUE(density)
  if (!is.null(aE)) p$aE <- as.numeric(aE)
  if (!is.null(xdist)) p$xdist <- isTRUE(xdist)
  if (!is.null(xdist_uniform)) p$xdist_uniform <- isTRUE(xdist_uniform)
  if (!is.null(xl1)) p$xl1 <- as.numeric(xl1)
  if (!is.null(tdist)) {
    p$tdist <- rep_len(as.logical(tdist), state$q_t)
  }
  if (!is.null(tdist_uniform)) {
    p$tdist_uniform <- rep_len(as.logical(tdist_uniform), state$q_t)
  }
  if (!is.null(tl1)) {
    p$tl1 <- rep_len(as.numeric(tl1), state$q_t)
  }
  if (!is.null(normal)) p$normal <- isTRUE(normal)
  if (!is.null(mu)) p$mu <- as.numeric(mu)
  if (!is.null(sigma)) p$sigma <- as.numeric(sigma)
  p
}

#' @noRd
#' @export
dynamics_description_length <- function(state, sbm = NULL, latent_edges = TRUE,
                                        density = FALSE, aE = 1,
                                        xdist = NULL, tdist = NULL,
                                        xdist_uniform = FALSE, tdist_uniform = FALSE,
                                        xl1 = NULL, tl1 = NULL,
                                        normal = FALSE, mu = 0, sigma = 1,
                                        breakdown = FALSE) {
  .check_dynamics_state(state)
  opts <- .dyn_merge_prior(
    state, sbm = sbm, latent_edges = latent_edges, density = density, aE = aE,
    xdist = xdist, tdist = tdist, xdist_uniform = xdist_uniform,
    tdist_uniform = tdist_uniform, xl1 = xl1, tl1 = tl1,
    normal = normal, mu = mu, sigma = sigma
  )
  res <- cpp_dynamics_entropy(
    state$model_code, state$n_vertices, state$directed, state$self_loops,
    state$s_list, state$f_mat, state$has_zero, state$exposed,
    state$positive, state$pslack, state$x_range[1L], state$x_range[2L],
    state$tmin, state$tmax, state$edges, state$x, state$theta, state$bs,
    state$degree_corrected, state$dist, opts, state$q_max
  )
  if (isTRUE(breakdown)) return(res)
  as.numeric(res$total)
}

#' @noRd
#' @export
reconstruct_dynamics <- function(state, niter = 25L, nsamples = 0L, burnin = 50L,
                                 thin = 2L, beta = Inf, k_cand = 20L,
                                 w_step = 0.25, t_step = 0.25, tol = 1e-4,
                                 threshold = 0.5, sbm = NULL,
                                 xdist = NULL, tdist = NULL,
                                 xl1 = NULL, tl1 = NULL,
                                 normal = FALSE, mu = 0, sigma = 1,
                                 refresh_every = 3L, merge_split = TRUE) {
  .check_dynamics_state(state)
  niter <- as.integer(niter)
  nsamples <- as.integer(nsamples)
  burnin <- as.integer(burnin)
  thin <- as.integer(thin)
  k_cand <- as.integer(k_cand)
  refresh_every <- as.integer(refresh_every)
  merge_split <- isTRUE(merge_split)
  if (is.na(niter) || niter < 0L) stop("niter must be >= 0", call. = FALSE)
  if (is.na(nsamples) || nsamples < 0L) stop("nsamples must be >= 0", call. = FALSE)
  if (is.na(refresh_every) || refresh_every < 0L) stop("refresh_every must be >= 0", call. = FALSE)
  if (niter == 0L && nsamples == 0L) stop("At least one of niter or nsamples must be > 0", call. = FALSE)

  opts <- .dyn_merge_prior(
    state, sbm = sbm, xdist = xdist, tdist = tdist,
    xl1 = xl1, tl1 = tl1, normal = normal, mu = mu, sigma = sigma
  )

  res <- cpp_dynamics_reconstruct(
    state$model_code, state$n_vertices, state$directed, state$self_loops,
    state$s_list, state$f_mat, state$has_zero, state$exposed,
    state$positive, state$pslack, state$x_range[1L], state$x_range[2L],
    state$tmin, state$tmax, state$edges, state$x, state$theta, state$bs,
    state$degree_corrected, state$dist, opts,
    niter, nsamples, burnin, thin, as.numeric(beta),
    k_cand, as.numeric(w_step), as.numeric(t_step), as.numeric(tol), state$q_max,
    refresh_every, merge_split
  )

  edges_df <- data.frame(
    from = res$from,
    to = res$to,
    prob = res$prob,
    weight = res$weight,
    weight_sd = res$weight_sd,
    stringsAsFactors = FALSE
  )
  if (nrow(edges_df) > 0L) {
    ord <- order(-edges_df$prob, -abs(edges_df$weight), edges_df$from, edges_df$to)
    edges_df <- edges_df[ord, , drop = FALSE]
    rownames(edges_df) <- NULL
  }

  g_rec <- graph_create(state$n_vertices, directed = state$directed)
  keep <- edges_df$prob >= threshold
  if (any(keep)) {
    sub_df <- edges_df[keep, , drop = FALSE]
    emat <- cbind(as.integer(sub_df$from), as.integer(sub_df$to))
    graph_add_edges(g_rec, emat, properties = list(
      weight = sub_df$weight,
      prob = sub_df$prob
    ))
  }

  new_state <- state
  new_state$edges <- res$final_edges
  new_state$x <- res$final_x
  new_state$theta <- res$final_theta
  new_state$bs <- res$final_bs
  new_state$default_prior <- opts

  structure(
    list(
      model = state$model,
      edges = edges_df,
      graph = g_rec,
      theta = res$theta_mean,
      theta_sd = res$theta_sd,
      bs_samples = if (length(res$bs_samples) > 0L) res$bs_samples else NULL,
      state = new_state,
      prior = opts,
      entropies = res$entropies,
      entropy = res$entropy,
      nll = res$nll,
      dS = res$dS,
      nattempts = res$nattempts,
      nmoves = res$nmoves
    ),
    class = "rgraphtool_dynamics_reconstruction"
  )
}

#' @noRd
#' @export
reconstruct_l1_path <- function(state, xl1_values = c(0.25, 0.5, 1, 2, 5, 10),
                                tl1 = 0.1, sbm = FALSE, niter = 20L,
                                truth = NULL, k_cand = 20L, tol = 1e-4) {
  .check_dynamics_state(state)
  if (!is.numeric(xl1_values) || length(xl1_values) < 1L || anyNA(xl1_values) || any(xl1_values < 0)) {
    stop("xl1_values must be a non-empty vector of non-negative numbers", call. = FALSE)
  }
  xl1_sorted <- sort(unique(as.numeric(xl1_values)), decreasing = TRUE)
  K <- length(xl1_sorted)
  fits <- vector("list", K)
  cur_state <- state

  # Prepare ground truth edge set if provided
  truth_keys <- NULL
  all_pairs <- NULL
  ek <- function(u, v, dir) {
    if (!dir) paste0(pmin(u, v), ":", pmax(u, v)) else paste0(u, ":", v)
  }
  if (!is.null(truth)) {
    if (inherits(truth, "rgraphtool_graph")) {
      te <- graph_edges(truth)
    } else {
      te <- .uncertain_validate_edges_mat(truth, state$n_vertices, "truth")
    }
    truth_keys <- unique(ek(te[, 1L], te[, 2L], state$directed))
    N <- state$n_vertices
    if (state$directed) {
      g_all <- expand.grid(u = seq_len(N), v = seq_len(N))
      if (!state$self_loops) g_all <- g_all[g_all$u != g_all$v, , drop = FALSE]
      all_pairs <- cbind(g_all$u, g_all$v)
    } else {
      all_pairs <- t(utils::combn(N, 2L))
      if (state$self_loops) all_pairs <- rbind(all_pairs, cbind(seq_len(N), seq_len(N)))
    }
  }

  sum_df <- data.frame(
    xl1 = xl1_sorted,
    n_edges = integer(K),
    nll = numeric(K),
    entropy = numeric(K),
    bic = numeric(K),
    l1_norm = numeric(K),
    stringsAsFactors = FALSE
  )
  if (!is.null(truth_keys)) {
    sum_df$precision <- NA_real_
    sum_df$recall <- NA_real_
    sum_df$f1 <- NA_real_
    sum_df$roc_auc <- NA_real_
  }

  n_obs <- state$n_obs
  for (i in seq_len(K)) {
    lam <- xl1_sorted[i]
    rec <- reconstruct_dynamics(
      cur_state, niter = niter, nsamples = 0L, beta = Inf,
      k_cand = k_cand, tol = tol,
      sbm = sbm, xdist = FALSE, tdist = FALSE, xl1 = lam, tl1 = tl1
    )
    fits[[i]] <- rec
    cur_state <- rec$state
    ne <- nrow(rec$edges)
    sum_df$n_edges[i] <- ne
    sum_df$nll[i] <- rec$nll
    sum_df$entropy[i] <- rec$entropy
    sum_df$bic[i] <- 2 * rec$nll + (ne + state$n_vertices * state$q_t) * log(max(2, n_obs))
    sum_df$l1_norm[i] <- sum(abs(rec$edges$weight))

    if (!is.null(truth_keys)) {
      pk <- ek(all_pairs[, 1L], all_pairs[, 2L], state$directed)
      y_true <- as.integer(pk %in% truth_keys)
      rk <- ek(rec$edges$from, rec$edges$to, state$directed)
      idx <- match(pk, rk)
      scores <- ifelse(is.na(idx), 0, abs(rec$edges$weight[idx]))
      pred_bin <- as.integer(!is.na(idx))
      tp <- sum(pred_bin == 1L & y_true == 1L)
      fp <- sum(pred_bin == 1L & y_true == 0L)
      fn <- sum(pred_bin == 0L & y_true == 1L)
      prec <- if (tp + fp > 0L) tp / (tp + fp) else NA_real_
      recl <- if (tp + fn > 0L) tp / (tp + fn) else NA_real_
      sum_df$precision[i] <- prec
      sum_df$recall[i] <- recl
      sum_df$f1[i] <- if (!is.na(prec) && !is.na(recl) && (prec + recl) > 0) {
        2 * prec * recl / (prec + recl)
      } else NA_real_
      if (sum(y_true == 1L) > 0L && sum(y_true == 0L) > 0L) {
        sum_df$roc_auc[i] <- evaluate_edge_prediction(scores, y_true)$roc_auc
      }
    }
  }

  best_idx <- which.min(sum_df$bic)
  structure(
    list(
      summary = sum_df,
      fits = fits,
      best_xl1 = sum_df$xl1[best_idx],
      best_fit = fits[[best_idx]]
    ),
    class = "rgraphtool_l1_path"
  )
}

#' @noRd
#' @export
dynamics_precision_matrix <- function(state) {
  if (inherits(state, "rgraphtool_dynamics_reconstruction")) state <- state$state
  .check_dynamics_state(state)
  if (!identical(state$model, "pseudo_normal")) {
    stop("dynamics_precision_matrix() requires a pseudo_normal_block_state()", call. = FALSE)
  }
  N <- state$n_vertices
  W <- matrix(0, nrow = N, ncol = N)
  if (nrow(state$edges) > 0L) {
    for (i in seq_len(nrow(state$edges))) {
      u <- state$edges[i, 1L]; v <- state$edges[i, 2L]; w <- state$x[i]
      W[u, v] <- w
      W[v, u] <- w
    }
  }
  th0 <- state$theta[, 1L]
  if (isTRUE(state$positive)) {
    k_abs <- rowSums(abs(W))
    idx <- k_abs > 0
    th0[idx] <- pmin(-0.5 * log(k_abs[idx]) - state$pslack, th0[idx])
  }
  diag(W) <- exp(-2 * th0)
  d_inv <- 1 / sqrt(diag(W))
  pcor <- -W * outer(d_inv, d_inv)
  diag(pcor) <- 1
  list(precision = W, partial_correlation = pcor, sigma = exp(th0), mu = state$theta[, 2L])
}

#' @noRd
#' @export
dynamics_transmission_probs <- function(state) {
  if (inherits(state, "rgraphtool_dynamics_reconstruction")) state <- state$state
  .check_dynamics_state(state)
  if (!identical(state$model, "epidemics")) {
    stop("dynamics_transmission_probs() requires an epidemics_block_state()", call. = FALSE)
  }
  edf <- if (nrow(state$edges) > 0L) {
    data.frame(
      from = state$edges[, 1L],
      to = state$edges[, 2L],
      x = state$x,
      transmission_prob = -expm1(state$x),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(from = integer(0), to = integer(0), x = numeric(0),
               transmission_prob = numeric(0))
  }
  list(edges = edf, spontaneous_prob = exp(state$theta[, 1L]))
}

# ---- Native C++ forward dynamical simulators ----

.sim_resolve_edge_param <- function(graph, val, default_val, name) {
  m <- graph_num_edges(graph)
  if (is.null(val)) return(rep(as.numeric(default_val), m))
  if (is.character(val) && length(val) == 1L) return(as.numeric(graph_edge_property(graph, val)))
  if (is.numeric(val) && length(val) == 1L) return(rep(as.numeric(val), m))
  if (is.numeric(val) && length(val) == m && !anyNA(val)) return(as.numeric(val))
  stop(sprintf("%s must be NULL, scalar, edge property name, or numeric vector of length %d",
               name, m), call. = FALSE)
}

.sim_resolve_node_param <- function(N, val, default_val, name) {
  if (is.null(val)) return(rep(as.numeric(default_val), N))
  if (is.numeric(val) && length(val) == 1L && !is.na(val)) return(rep(as.numeric(val), N))
  if (is.numeric(val) && length(val) == N && !anyNA(val)) return(as.numeric(val))
  stop(sprintf("%s must be NULL, scalar, or numeric vector of length %d", name, N), call. = FALSE)
}

#' @noRd
#' @export
simulate_ising <- function(graph, weights = 0.5, theta = 0, n_samples = 200L,
                           burnin = 100L, thin = 1L,
                           mode = c("equilibrium", "glauber"),
                           has_zero = FALSE, init_state = NULL) {
  ptr <- graph_pointer(graph)
  N <- graph_num_vertices(graph)
  mode <- match.arg(mode)
  w <- .sim_resolve_edge_param(graph, weights, 0.5, "weights")
  th <- .sim_resolve_node_param(N, theta, 0, "theta")
  if (is.null(init_state)) {
    s0 <- sample(if (isTRUE(has_zero)) c(-1L, 0L, 1L) else c(-1L, 1L), N, replace = TRUE)
  } else {
    s0 <- as.integer(init_state)
    if (length(s0) != N || anyNA(s0)) stop("init_state must have length N", call. = FALSE)
  }
  cpp_sim_ising(ptr, w, th, as.integer(n_samples), as.integer(burnin),
                as.integer(thin), identical(mode, "equilibrium"),
                isTRUE(has_zero), s0)
}

#' @noRd
#' @export
simulate_epidemics <- function(graph, transmission = 0.3, spontaneous = 0.005,
                               recovery = 0.2, resusceptible = 0.0,
                               model = c("SI", "SIS", "SIR", "SIRS"),
                               n_steps = 50L, n_cascades = 1L,
                               exposed = FALSE, expose_rate = 0.5,
                               init_infected = NULL) {
  ptr <- graph_pointer(graph)
  N <- graph_num_vertices(graph)
  model <- match.arg(model)
  t_prob <- .sim_resolve_edge_param(graph, transmission, 0.3, "transmission")
  if (any(t_prob < 0 | t_prob >= 1)) stop("transmission probabilities must lie in [0, 1)", call. = FALSE)
  log1m_t <- log1p(-t_prob)
  sp <- .sim_resolve_node_param(N, spontaneous, 0.005, "spontaneous")
  if (any(sp < 0 | sp >= 1)) stop("spontaneous infection probabilities must lie in [0, 1)", call. = FALSE)

  gam <- if (model == "SI") rep(0, N) else .sim_resolve_node_param(N, recovery, 0.2, "recovery")
  mu_r <- if (model == "SIS") {
    rep(-1, N)
  } else if (model == "SIRS") {
    .sim_resolve_node_param(N, resusceptible, 0.1, "resusceptible")
  } else {
    rep(0, N)
  }

  n_cascades <- as.integer(n_cascades)
  n_steps <- as.integer(n_steps)
  run_one <- function() {
    s0 <- integer(N)
    if (is.null(init_infected)) {
      s0[sample.int(N, 1L)] <- 1L
    } else if (length(init_infected) == N && all(init_infected %in% c(-1L, 0L, 1L, 2L))) {
      s0 <- as.integer(init_infected)
    } else {
      s0[as.integer(init_infected)] <- 1L
    }
    cpp_sim_epidemics(ptr, log1m_t, sp, gam, mu_r, n_steps,
                      isTRUE(exposed), as.numeric(expose_rate), s0)
  }
  if (n_cascades == 1L) return(run_one())
  replicate(n_cascades, run_one(), simplify = FALSE)
}

#' @noRd
#' @export
simulate_normal_dynamics <- function(graph, weights = 0.1, sigma = 1, mu = 0,
                                     n_samples = 200L, burnin = 100L, thin = 1L,
                                     mode = c("equilibrium", "glauber", "linear")) {
  ptr <- graph_pointer(graph)
  N <- graph_num_vertices(graph)
  mode <- match.arg(mode)
  w <- .sim_resolve_edge_param(graph, weights, 0.1, "weights")
  sig <- .sim_resolve_node_param(N, sigma, 1, "sigma")
  if (any(sig <= 0)) stop("sigma must be positive", call. = FALSE)
  mu_v <- .sim_resolve_node_param(N, mu, 0, "mu")
  mode_code <- match(mode, c("equilibrium", "glauber", "linear")) - 1L
  cpp_sim_normal_dyn(ptr, w, sig, mu_v, as.integer(n_samples),
                     as.integer(burnin), as.integer(thin), as.integer(mode_code))
}

#' @noRd
#' @export
simulate_potts <- function(graph, weights = 0.5, f = diag(3), theta = NULL,
                           n_samples = 200L, burnin = 100L, thin = 1L,
                           mode = c("equilibrium", "glauber"), init_state = NULL) {
  ptr <- graph_pointer(graph)
  N <- graph_num_vertices(graph)
  mode <- match.arg(mode)
  f_mat <- as.matrix(f)
  storage.mode(f_mat) <- "double"
  q <- nrow(f_mat)
  w <- .sim_resolve_edge_param(graph, weights, 0.5, "weights")
  th <- .dyn_resolve_theta(theta, N, q)
  s0 <- if (is.null(init_state)) sample(0:(q - 1L), N, replace = TRUE) else as.integer(init_state)
  cpp_sim_potts(ptr, w, f_mat, th, as.integer(n_samples), as.integer(burnin),
                as.integer(thin), identical(mode, "equilibrium"), s0)
}

# ---- Print methods ----

#' @noRd
#' @export
print.rgraphtool_dynamics_block_state <- function(x, ...) {
  cat(sprintf("<rgraphtool_dynamics_block_state (%s, %s, %s)>\n",
              x$model,
              if (x$directed) "directed" else "undirected",
              if (x$nested) "nested" else "flat"))
  cat(sprintf("  Vertices: %d | Observations: %d | Edges: %d | Node params/vertex: %d\n",
              x$n_vertices, x$n_obs, nrow(x$edges), x$q_t))
  cat(sprintf("  Prior: sbm=%s, xdist=%s, xl1=%.3g, tdist=%s\n",
              x$default_prior$sbm, x$default_prior$xdist, x$default_prior$xl1,
              any(x$default_prior$tdist)))
  invisible(x)
}

#' @noRd
#' @export
print.rgraphtool_dynamics_reconstruction <- function(x, ...) {
  cat(sprintf("<rgraphtool_dynamics_reconstruction (%s)>\n", x$model))
  cat(sprintf("  Inferred edges: %d | Final DL: %.4f (NLL: %.4f)\n",
              nrow(x$edges), x$entropy, x$nll))
  cat(sprintf("  Moves: %.0f / %.0f | Prior: sbm=%s, xdist=%s, xl1=%.3g\n",
              x$nmoves, x$nattempts, x$prior$sbm, x$prior$xdist, x$prior$xl1))
  invisible(x)
}

#' @noRd
#' @export
print.rgraphtool_l1_path <- function(x, ...) {
  cat(sprintf("<rgraphtool_l1_path (%d lambda values, best BIC at xl1 = %.4g)>\n",
              nrow(x$summary), x$best_xl1))
  print(x$summary, row.names = FALSE)
  invisible(x)
}
