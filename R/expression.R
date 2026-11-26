# Gene-expression preprocessing, exact fixed-effects QR projection,
# design-aware gene filtering, and SBM module eigengene helpers.

.expr_check_matrix <- function(expr) {
  if (is.data.frame(expr)) expr <- as.matrix(expr)
  if (!is.matrix(expr) || !is.numeric(expr) || nrow(expr) < 1L || ncol(expr) < 2L ||
      anyNA(expr) || any(!is.finite(expr))) {
    stop("expr must be a finite numeric matrix with at least 1 row (gene) and 2 columns (samples)",
         call. = FALSE)
  }
  storage.mode(expr) <- "double"
  expr
}

.expr_design_matrix <- function(design, col_data, n_samples) {
  if (inherits(design, "formula")) {
    if (is.null(col_data) || !is.data.frame(col_data)) {
      stop("col_data must be a data.frame when design is a formula", call. = FALSE)
    }
    if (nrow(col_data) != n_samples) {
      stop(sprintf("col_data has %d rows, expected %d matching ncol(expr)",
                   nrow(col_data), n_samples), call. = FALSE)
    }
    X <- stats::model.matrix(design, data = col_data)
  } else {
    X <- as.matrix(design)
    if (!is.numeric(X) || nrow(X) != n_samples || anyNA(X) || any(!is.finite(X))) {
      stop(sprintf("design matrix must be finite numeric with %d rows matching ncol(expr)",
                   n_samples), call. = FALSE)
    }
    storage.mode(X) <- "double"
  }
  X
}

#' Exact QR / REML fixed-effects projection for Gaussian network reconstruction
#'
#' Projects an `N x M` gene-expression matrix onto an orthonormal basis `Q` of
#' the orthogonal complement of the `M x p` design matrix `X`, producing an
#' `N x (M - r)` matrix `Z = expr %*% Q` whose columns are exact i.i.d. draws
#' with exact residual degrees of freedom `df_resid = M - r`.
#'
#' @param expr Numeric `N x M` matrix (genes x samples).
#' @param design A model formula (evaluated against `col_data`) or `M x p`
#'   numeric design matrix.
#' @param col_data Optional `data.frame` of sample metadata when `design` is a formula.
#' @param transform Residual scaling/transformation: `"zscore"` (default,
#'   scales each gene by residual SD), `"nonparanormal"` (rank-Gaussian copula),
#'   or `"none"` (raw projected residuals).
#' @return An S3 object of class `"rgraphtool_expression_projection"` containing
#'   `residuals`, `raw_residuals`, `beta`, `resid_var`, `r2_fixed`, `df_resid`,
#'   `Q`, `design_matrix`, and `transform`.
#' @export
project_fixed_effects <- function(expr, design, col_data = NULL,
                                  transform = c("zscore", "nonparanormal", "none")) {
  transform <- match.arg(transform)
  expr <- .expr_check_matrix(expr)
  M <- ncol(expr)
  N <- nrow(expr)

  X <- .expr_design_matrix(design, col_data, M)
  qr_x <- qr(X)
  r <- qr_x$rank
  if (r >= M) {
    stop(sprintf("Design rank (%d) must be strictly less than ncol(expr) (%d)", r, M),
         call. = FALSE)
  }

  Q_full <- qr.Q(qr_x, complete = TRUE)
  Q <- Q_full[, (r + 1L):M, drop = FALSE]
  df_resid <- M - r

  Z_raw <- expr %*% Q
  beta <- t(qr.coef(qr_x, t(expr)))
  if (!is.matrix(beta)) {
    beta <- matrix(beta, nrow = N, ncol = ncol(X))
  }
  rownames(beta) <- rownames(expr)
  colnames(beta) <- colnames(X)

  rss <- rowSums(Z_raw^2)
  resid_var <- rss / df_resid
  tot_ss <- rowSums((expr - rowMeans(expr))^2)
  r2_fixed <- ifelse(tot_ss > 0, pmax(0, 1 - rss / tot_ss), 0)

  if (transform == "zscore") {
    # Ensure exact zero mean if X lacked an intercept, and unit sample variance
    row_means <- rowMeans(Z_raw)
    Z_c <- Z_raw - row_means
    row_sd <- sqrt(pmax(rowSums(Z_c^2) / max(1L, df_resid - 1L), 1e-12))
    Z <- Z_c / row_sd
  } else if (transform == "nonparanormal") {
    Z <- t(apply(Z_raw, 1L, function(row) {
      rk <- rank(row, ties.method = "average")
      z <- stats::qnorm(rk / (length(row) + 1))
      s <- stats::sd(z)
      if (is.na(s) || s <= 1e-12) rep(0, length(z)) else (z - mean(z)) / s
    }))
    if (N == 1L) Z <- matrix(Z, nrow = 1L)
  } else {
    Z <- Z_raw
  }
  rownames(Z) <- rownames(expr)
  rownames(Z_raw) <- rownames(expr)

  structure(
    list(
      residuals = Z,
      raw_residuals = Z_raw,
      beta = beta,
      resid_var = resid_var,
      r2_fixed = r2_fixed,
      df_resid = df_resid,
      rank = r,
      Q = Q,
      design_matrix = X,
      transform = transform
    ),
    class = "rgraphtool_expression_projection"
  )
}

#' @noRd
#' @export
print.rgraphtool_expression_projection <- function(x, ...) {
  cat(sprintf(
    "<expression projection: %d genes, %d samples -> %d orthogonal df (rank %d removed, transform = '%s')>\n",
    nrow(x$residuals), nrow(x$design_matrix), x$df_resid, x$rank, x$transform
  ))
  cat(sprintf(
    "  median R^2 explained by fixed effects: %.3f\n",
    stats::median(x$r2_fixed)
  ))
  invisible(x)
}

#' Design-aware gene filtering via residual variance and nested-model ANOVA
#'
#' @param expr Numeric `N x M` expression matrix.
#' @param col_data Sample metadata `data.frame` (required if `nuisance` or
#'   `signal` is a formula).
#' @param nuisance Optional nuisance formula or matrix (defaults to intercept-only).
#' @param signal Optional signal formula or matrix tested beyond `nuisance` via
#'   an exact vectorized F-test.
#' @param min_var Minimum residual variance after projecting out `nuisance`.
#' @param max_padj Maximum Benjamini-Hochberg adjusted p-value (only applied
#'   when `signal` is provided).
#' @param top_n Optional cap on retained genes (ranked by F-statistic if
#'   `signal` is given, otherwise by residual variance).
#' @return List with filtered `expr`, `stats` data.frame, and logical `keep`.
#' @export
filter_expression_genes <- function(expr, col_data = NULL,
                                    nuisance = NULL, signal = NULL,
                                    min_var = 0, max_padj = 1, top_n = NULL) {
  expr <- .expr_check_matrix(expr)
  N <- nrow(expr)
  M <- ncol(expr)

  X_nuis <- if (is.null(nuisance)) {
    matrix(1, nrow = M, ncol = 1L)
  } else {
    .expr_design_matrix(nuisance, col_data, M)
  }
  qr_nuis <- qr(X_nuis)
  r_nuis <- qr_nuis$rank
  if (r_nuis >= M) {
    stop("Nuisance design leaves zero residual degrees of freedom", call. = FALSE)
  }
  Q_nuis <- qr.Q(qr_nuis, complete = TRUE)[, (r_nuis + 1L):M, drop = FALSE]
  Z_nuis <- expr %*% Q_nuis
  rss_nuis <- rowSums(Z_nuis^2)
  df_nuis <- M - r_nuis
  var_nuis <- rss_nuis / df_nuis

  tot_ss <- rowSums((expr - rowMeans(expr))^2)
  r2_nuis <- ifelse(tot_ss > 0, pmax(0, 1 - rss_nuis / tot_ss), 0)

  gene_ids <- if (!is.null(rownames(expr))) rownames(expr) else paste0("G", seq_len(N))
  stats_df <- data.frame(
    gene = gene_ids,
    var_resid = var_nuis,
    r2_nuis = r2_nuis,
    stringsAsFactors = FALSE
  )

  has_signal <- !is.null(signal)
  if (has_signal) {
    X_sig <- .expr_design_matrix(signal, col_data, M)
    X_full <- cbind(X_nuis, X_sig)
    qr_full <- qr(X_full)
    r_full <- qr_full$rank
    df1 <- r_full - r_nuis
    df2 <- M - r_full
    if (df1 <= 0L || df2 <= 0L) {
      stop(sprintf("Nested model test requires df1 > 0 and df2 > 0 (got df1 = %d, df2 = %d)",
                   df1, df2), call. = FALSE)
    }
    Q_full <- qr.Q(qr_full, complete = TRUE)[, (r_full + 1L):M, drop = FALSE]
    rss_full <- rowSums((expr %*% Q_full)^2)
    f_stat <- ((pmax(0, rss_nuis - rss_full)) / df1) / pmax(rss_full / df2, 1e-15)
    pval <- stats::pf(f_stat, df1, df2, lower.tail = FALSE)
    padj <- stats::p.adjust(pval, method = "BH")
    stats_df$F_stat <- f_stat
    stats_df$pvalue <- pval
    stats_df$padj <- padj
  }

  keep <- stats_df$var_resid >= min_var
  if (has_signal) {
    keep <- keep & !is.na(stats_df$padj) & (stats_df$padj <= max_padj)
  }

  if (!is.null(top_n)) {
    top_n <- as.integer(top_n)
    cand <- which(keep)
    if (length(cand) > top_n) {
      score <- if (has_signal) stats_df$F_stat[cand] else stats_df$var_resid[cand]
      chosen <- cand[order(score, decreasing = TRUE)[seq_len(top_n)]]
      keep <- rep(FALSE, N)
      keep[chosen] <- TRUE
    }
  }
  stats_df$keep <- keep

  list(
    expr = expr[keep, , drop = FALSE],
    stats = stats_df,
    keep = keep
  )
}

#' Compute sign-aligned module eigengenes (PC1) for SBM blocks
#'
#' @param expr Numeric `N x M` expression matrix.
#' @param blocks Integer block vector of length `N`, or an `rgraphtool` flat/nested
#'   block state.
#' @param level Hierarchy level (`0L` = finest) when `blocks` is a nested state.
#' @param min_size Minimum number of genes required to compute a module eigengene.
#' @return List with `eigengenes` (`M x B_kept` matrix), `pve` (variance
#'   explained by PC1 per module), and `block_sizes`.
#' @export
module_eigengenes <- function(expr, blocks, level = 0L, min_size = 2L) {
  expr <- .expr_check_matrix(expr)
  if (is_nested_state(blocks)) {
    b_vec <- nested_project(blocks, level = level)
  } else if (inherits(blocks, "rgraphtool_block_state") ||
             inherits(blocks, "rgraphtool_weighted_block_state") ||
             inherits(blocks, "rgraphtool_layered_block_state") ||
             inherits(blocks, "rgraphtool_pp_block_state") ||
             inherits(blocks, "rgraphtool_ranked_block_state")) {
    b_vec <- block_membership(blocks)
  } else {
    b_vec <- as.integer(blocks)
  }
  if (length(b_vec) != nrow(expr)) {
    stop("blocks length must match nrow(expr)", call. = FALSE)
  }

  # Standardize each gene across samples
  z_expr <- t(scale(t(expr)))
  z_expr[!is.finite(z_expr)] <- 0

  tab <- table(b_vec)
  kept_blocks <- as.integer(names(tab)[tab >= min_size])
  if (length(kept_blocks) == 0L) {
    return(list(
      eigengenes = matrix(numeric(0), nrow = ncol(expr), ncol = 0L,
                          dimnames = list(colnames(expr), character(0))),
      pve = numeric(0),
      block_sizes = integer(0)
    ))
  }

  me_mat <- matrix(0, nrow = ncol(expr), ncol = length(kept_blocks))
  pve <- numeric(length(kept_blocks))
  sizes <- integer(length(kept_blocks))

  for (i in seq_along(kept_blocks)) {
    b <- kept_blocks[i]
    idx <- which(b_vec == b)
    sizes[i] <- length(idx)
    sub_mat <- z_expr[idx, , drop = FALSE]
    sv <- svd(t(sub_mat), nu = 1L, nv = 0L)
    pc1 <- as.numeric(sv$u[, 1L]) * sv$d[1L]
    avg_expr <- colMeans(sub_mat)
    if (stats::cor(pc1, avg_expr) < 0) pc1 <- -pc1
    s_pc <- stats::sd(pc1)
    if (!is.na(s_pc) && s_pc > 0) pc1 <- (pc1 - mean(pc1)) / s_pc
    me_mat[, i] <- pc1
    tot_var <- sum(sv$d^2)
    pve[i] <- if (tot_var > 0) (sv$d[1L]^2) / tot_var else 0
  }

  col_lbls <- paste0("ME", kept_blocks)
  dimnames(me_mat) <- list(colnames(expr), col_lbls)
  names(pve) <- col_lbls
  names(sizes) <- col_lbls

  list(
    eigengenes = me_mat,
    pve = pve,
    block_sizes = sizes
  )
}

#' Analytical Ledoit-Wolf / Schafer-Strimmer shrinkage covariance estimation
#'
#' Computes an optimal linear shrinkage covariance matrix
#' \eqn{(1 - \lambda) S + \lambda T} for high-dimensional expression data
#' (`N >> M`) with zero external package dependencies, suitable for passing
#' directly to `pseudo_normal_block_state(cov = ..., n_eff = ...)`.
#'
#' @param expr Numeric `N x M` matrix (genes x samples).
#' @param target Shrinkage target: `"diagonal"` (unequal gene variances,
#'   zero covariances) or `"identity"` (mean-variance scaled identity).
#' @return List with `cov` (`N x N` symmetric SPD matrix), `lambda` (shrinkage
#'   intensity in `[0, 1]`), and `n_obs` (`M`).
#' @export
shrinkage_covariance <- function(expr, target = c("diagonal", "identity")) {
  target <- match.arg(target)
  expr <- .expr_check_matrix(expr)
  N <- nrow(expr)
  M <- ncol(expr)
  if (M < 3L) {
    stop("shrinkage_covariance() requires at least 3 samples (columns)", call. = FALSE)
  }

  xc <- expr - rowMeans(expr)
  s_mat <- tcrossprod(xc) / (M - 1)
  dimnames(s_mat) <- list(rownames(expr), rownames(expr))

  if (N == 1L) {
    return(list(cov = s_mat, lambda = 0, n_obs = M))
  }

  # Unbiased variance of sample covariances (Schafer & Strimmer 2005)
  z_bar <- ((M - 1) / M) * s_mat
  w2_sum <- tcrossprod(xc * xc)
  diff_w <- w2_sum - M * (z_bar * z_bar)
  diff_w[diff_w < 0] <- 0
  var_s <- (M / ((M - 1)^3)) * diff_w

  if (target == "diagonal") {
    t_mat <- diag(diag(s_mat), nrow = N, ncol = N)
    sum_var <- sum(var_s) - sum(diag(var_s))
    sum_sq <- sum(s_mat * s_mat) - sum(diag(s_mat)^2)
  } else {
    v_mean <- mean(diag(s_mat))
    t_mat <- diag(v_mean, nrow = N, ncol = N)
    sum_var <- sum(var_s)
    diff_mat <- s_mat - t_mat
    sum_sq <- sum(diff_mat * diff_mat)
  }

  lam <- if (sum_sq > 0) min(1, max(0, sum_var / sum_sq)) else 0
  cov_sh <- (1 - lam) * s_mat + lam * t_mat
  cov_sh <- 0.5 * (cov_sh + t(cov_sh))
  dimnames(cov_sh) <- list(rownames(expr), rownames(expr))

  list(cov = cov_sh, lambda = lam, n_obs = M)
}
