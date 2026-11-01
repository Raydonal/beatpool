#' Extract PCA bases from a list of beat matrices
#'
#' @param X_list A list of beat matrices, one per record.
#' @param q Number of components to keep.
#' @param ... Passed to [pca_basis()].
#' @return A list of \code{bp_basis} objects.
#' @export
extract_bases <- function(X_list, q = 10L, ...) {
  if (!is.list(X_list)) stop("X_list must be a list of matrices", call. = FALSE)
  lapply(X_list, function(X) pca_basis(X, q = q, ...))
}

.as_basis_matrix <- function(b) {
  M <- if (is.matrix(b)) b else if (inherits(b, "bp_basis")) b$T else
    if (is.list(b) && !is.null(b$T)) b$T else NULL
  if (is.null(M)) stop("cannot extract a basis matrix", call. = FALSE)
  .check_numeric_matrix(M, "basis", min_rows = 1L)
}

#' Pool patient-specific bases
#'
#' Combines a list of orthonormal bases into a single basis. Four pooling
#' operators are available: the arithmetic mean, the column-normalised mean,
#' the QR retraction and the orthogonal Procrustes mean. The Procrustes mean is
#' the minimiser of the sum of squared Frobenius distances to the inputs and is
#' the recommended default.
#'
#' @param bases A list of matrices or \code{bp_basis} objects.
#' @param method One of \code{"procrustes"}, \code{"bar"}, \code{"qr"} or
#'   \code{"mean"}.
#' @param align Align the bases before pooling. Alignment removes the sign,
#'   permutation and rotation indeterminacy of principal components.
#' @param ref Reference used by the one-step alignment: an integer index into
#'   \code{bases}, or a matrix. Ignored when \code{iterations > 1}.
#' @param iterations Number of alignment/pooling iterations. Values greater
#'   than one run generalised Procrustes analysis.
#' @param q Number of components to pool. Defaults to the smallest number of
#'   columns among the inputs.
#' @param data Optional list of beat matrices used to order the pooled columns
#'   by average explained variance.
#' @param weights Optional record weights.
#' @param tol Convergence tolerance for the iterative algorithm.
#' @param max_iter Maximum number of iterations.
#' @param robust Apply robust reweighting to downweight outlying records.
#' @param robust_c Tuning constant for the robust weights.
#' @return A list of class \code{bp_pool} with the pooled basis \code{T}, the
#'   method, the aligned bases, the per-record rotations, the average explained
#'   variances and convergence information.
#' @export
pool_bases <- function(bases, method = c("procrustes", "bar", "qr", "mean"),
                       align = TRUE, ref = 1L, iterations = 1L, q = NULL,
                       data = NULL, weights = NULL, tol = 1e-10,
                       max_iter = 50L, robust = FALSE, robust_c = 1.345) {
  method <- match.arg(method)
  if (inherits(bases, "bp_basis") || is.matrix(bases)) bases <- list(bases)
  if (!length(bases)) stop("no bases to pool", call. = FALSE)
  mats <- lapply(bases, .as_basis_matrix)
  if (is.null(q)) q <- min(vapply(mats, ncol, integer(1)))
  q <- as.integer(q)
  mats <- lapply(mats, function(M) M[, seq_len(min(q, ncol(M))), drop = FALSE])
  if (is.null(weights)) weights <- rep(1, length(mats))
  weights <- as.numeric(weights) / sum(weights)
  if (length(weights) != length(mats)) stop("weights have the wrong length",
                                            call. = FALSE)
  iterations <- max(1L, as.integer(iterations))

  mean_basis <- function(L) {
    A <- Reduce(`+`, Map(function(M, w) M * w, L, as.list(weights)))
    A
  }
  pool_mean <- function(A) {
    if (method == "mean") return(A)
    if (method == "bar") {
      nrm <- sqrt(colSums(A^2))
      if (any(nrm < .Machine$double.eps)) stop("degenerate pooled column",
                                               call. = FALSE)
      return(sweep(A, 2L, nrm, "/"))
    }
    if (method == "qr") return(qr.Q(qr(A)))
    s <- svd(A, nu = ncol(A), nv = ncol(A))
    s$u %*% t(s$v)
  }

  if (align) {
    converged <- NA
    it_used <- 1L
    if (is.character(ref) && length(ref) == 1L && ref == "medoid") {
      D <- matrix(0, length(mats), length(mats))
      for (i in seq_along(mats)) for (j in seq_along(mats)) {
        if (i < j) {
          D[i, j] <- D[j, i] <- subspace_distance(mats[[i]], mats[[j]])
        }
      }
      ref_idx <- which.min(rowSums(D))
      ref_mat <- mats[[ref_idx]]
    } else {
      ref_mat <- if (is.matrix(ref)) ref else mats[[as.integer(ref)]]
    }
    if (iterations == 1L) {
      aligned <- lapply(mats, function(M) align_basis(M, ref_mat)$Q)
      A <- mean_basis(aligned)
      Tmat <- sign_fix(pool_mean(A))
    } else {
      C <- ref_mat
      converged <- FALSE
      it_used <- 0L
      history <- vector("list", max_iter)
      changes <- numeric(max_iter)
      for (it in seq_len(max_iter)) {
        aligned <- lapply(mats, function(M) align_basis(M, C)$Q)
        Cnew <- mean_basis(aligned)
        d <- norm(Cnew - C, "F") / max(norm(C, "F"), .Machine$double.eps)
        C <- Cnew
        history[[it]] <- C
        changes[it] <- d
        it_used <- it
        if (d < tol) { converged <- TRUE; break }
      }
      used <- seq_len(it_used)
      best <- used[which.min(changes[used])]
      C <- history[[best]]
      it_used <- best
      converged <- changes[best] < sqrt(tol)
      Tmat <- sign_fix(pool_mean(C))
      aligned <- lapply(mats, function(M) align_basis(M, Tmat)$Q)
      attr(Tmat, "changes") <- changes[used]
    }
  } else {
    aligned <- mats
    Tmat <- sign_fix(pool_mean(mean_basis(mats)))
    converged <- NA
    it_used <- 1L
  }

  rotations <- lapply(mats, function(M) align_basis(M, Tmat)$O)
  effective_weights <- weights
  if (robust) {
    base_w <- weights
    for (r in seq_len(5L)) {
      d <- vapply(aligned, function(M) norm(M - Tmat, "F"), numeric(1))
      med <- stats::median(d)
      if (!is.finite(med) || med < .Machine$double.eps) break
      w <- base_w * pmin(1, robust_c * med / pmax(d, .Machine$double.eps))
      w <- w / sum(w)
      A <- Reduce(`+`, Map(function(M, ww) M * ww, aligned, as.list(w)))
      Tmat <- sign_fix(pool_mean(A))
      aligned <- lapply(mats, function(M) align_basis(M, Tmat)$Q)
      effective_weights <- w
    }
    rotations <- lapply(mats, function(M) align_basis(M, Tmat)$O)
  }
  values <- if (!is.null(data)) basis_variance(Tmat, data, weights) else NULL
  if (!is.null(values)) {
    o <- order(values, decreasing = TRUE)
    Tmat <- sign_fix(Tmat[, o, drop = FALSE])
    values <- values[o]
    aligned <- lapply(aligned, function(M) M[, o, drop = FALSE])
  }
  out <- list(
    T = Tmat, method = method, align = align, q = ncol(Tmat),
    aligned = aligned, rotations = rotations, values = values,
    weights = effective_weights, ref = ref, iterations = iterations,
    converged = converged, iterations_used = it_used,
    n_bases = length(mats)
  )
  class(out) <- c("bp_pool", "list")
  out
}

#' Generalised Procrustes analysis of bases
#'
#' Iteratively aligns and averages a list of bases until the pooled subspace
#' stabilises.
#'
#' @inheritParams pool_bases
#' @return A \code{bp_pool} object.
#' @export
gpa_bases <- function(bases, method = "procrustes", q = NULL, data = NULL,
                      weights = NULL, tol = 1e-10, max_iter = 50L) {
  pool_bases(bases, method = method, align = TRUE, iterations = max_iter,
             q = q, data = data, weights = weights, tol = tol,
             max_iter = max_iter)
}

#' Hierarchical Procrustes pooling
#'
#' Pools the leading \code{k}-dimensional subspaces separately for every
#' \code{k = 1, ..., q}. Each subspace is aligned to a medoid reference before
#' averaging, which avoids forcing a single compromise basis to serve all
#' dimensions.
#'
#' @param bases A list of matrices or \code{bp_basis} objects.
#' @param q Maximum number of components.
#' @param ref Reference rule: \code{"medoid"} (default) or an integer index.
#'   The medoid is computed once in the full \code{q}-dimensional subspace and
#'   the same reference record is used for every dimension.
#' @param iterations Number of GPA iterations for each subspace.
#' @param data Optional list of beat matrices used to order components.
#' @param robust Use robust reweighting.
#' @param ... Passed to [pool_bases()].
#' @return A list of class \code{bp_hpool} with one \code{bp_pool} object per
#'   dimension.
#' @export
hierarchical_pool <- function(bases, q = 10L, ref = "medoid", iterations = 1L,
                              data = NULL, robust = FALSE, ...) {
  if (inherits(bases, "bp_basis") || is.matrix(bases)) bases <- list(bases)
  mats <- lapply(bases, .as_basis_matrix)
  q <- min(as.integer(q), min(vapply(mats, ncol, integer(1))))
  if (is.character(ref) && length(ref) == 1L && ref == "medoid") {
    D <- matrix(0, length(mats), length(mats))
    for (i in seq_along(mats)) for (j in seq_along(mats)) {
      if (i < j) D[i, j] <- D[j, i] <- subspace_distance(mats[[i]], mats[[j]])
    }
    ref <- which.min(rowSums(D))
  }
  out <- vector("list", q)
  for (k in seq_len(q)) {
    b <- lapply(mats, function(M) M[, seq_len(k), drop = FALSE])
    d <- if (!is.null(data)) lapply(data, function(X) X) else NULL
    out[[k]] <- pool_bases(b, method = "procrustes", align = TRUE, ref = ref,
                           iterations = iterations, q = k, data = d,
                           robust = robust, ...)
  }
  names(out) <- paste0("k", seq_len(q))
  structure(out, class = c("bp_hpool", "list"), q = q)
}

#' @export
print.bp_hpool <- function(x, ...) {
  cat("beatpool hierarchical pooled bases\n")
  cat("  dimensions:", length(x), "| method: Procrustes\n")
  invisible(x)
}

#' PCA basis from an aggregated covariance matrix
#'
#' Computes the eigenvectors of the average record covariance matrix. This is
#' the basis obtained by a federated PCA that shares \code{p x p} covariance
#' matrices instead of \code{p x q} bases.
#'
#' @param X_list A list of beat matrices.
#' @param q Number of components.
#' @param weights Optional record weights.
#' @param center Center each record before computing its covariance.
#' @param return_cov Return the assembled covariance matrix as well.
#' @return A \code{bp_basis} object; with \code{return_cov = TRUE} a list with
#'   the basis \code{T} and the covariance \code{S}.
#' @export
pool_covariance <- function(X_list, q = 10L, weights = NULL, center = TRUE,
                            return_cov = FALSE) {
  if (!is.list(X_list)) X_list <- list(X_list)
  if (is.null(weights)) weights <- rep(1, length(X_list))
  weights <- weights / sum(weights)
  p <- ncol(X_list[[1]])
  S <- matrix(0, p, p)
  for (h in seq_along(X_list)) {
    X <- as.matrix(X_list[[h]])
    if (center) X <- sweep(X, 2L, colMeans(X), "-")
    S <- S + weights[h] * crossprod(X) / max(1, nrow(X) - 1L)
  }
  e <- eigen(S, symmetric = TRUE)
  q <- min(as.integer(q), p)
  out <- list(T = sign_fix(e$vectors[, seq_len(q), drop = FALSE]),
              values = e$values[seq_len(q)], mean = rep(0, p), q = q,
              n = sum(vapply(X_list, nrow, integer(1))), p = p,
              method = "covariance")
  if (return_cov) list(T = out$T, S = S) else out
}

#' Pool second moments from per-record spectral summaries
#'
#' Each record is summarised by its leading eigenpairs \eqn{(Q_h, \lambda_h)}
#' and the pooled covariance is assembled as
#' \deqn{\hat\Sigma = \sum_h w_h Q_h \Lambda_h Q_h'.}
#' The outer product is invariant to signs, permutations and rotations of
#' \eqn{Q_h}, so no alignment is needed. When \code{residual = TRUE} the
#' variance outside the retained subspace is spread over its orthogonal
#' complement, which keeps the trace of each record covariance intact.
#'
#' @param bases A list of \code{bp_basis} objects.
#' @param q Number of leading eigenpairs used in the per-record summaries.
#' @param q_out Number of components returned. It may exceed \code{q}.
#' @param weights Optional record weights.
#' @param residual Impute the residual variance over the discarded subspace.
#' @param return_cov Return the assembled covariance matrix as well.
#' @return A \code{bp_basis} object; with \code{return_cov = TRUE} a list with
#'   the basis \code{T} and the covariance \code{S}.
#' @export
pool_spectral <- function(bases, q = 10L, q_out = q, weights = NULL,
                          residual = FALSE, return_cov = FALSE) {
  if (inherits(bases, "bp_basis") || is.matrix(bases)) bases <- list(bases)
  if (is.null(weights)) weights <- rep(1, length(bases))
  weights <- weights / sum(weights)
  q <- as.integer(q)
  p <- nrow(.as_basis_matrix(bases[[1]]))
  S <- matrix(0, p, p)
  for (h in seq_along(bases)) {
    b <- bases[[h]]
    Q <- if (inherits(b, "bp_basis")) b$T else as.matrix(b)
    lam <- if (inherits(b, "bp_basis") && !is.null(b$values)) b$values else
      rep(1, ncol(Q))
    j <- seq_len(min(q, ncol(Q)))
    Q <- Q[, j, drop = FALSE]
    lam <- pmax(lam[j], 0)
    Sh <- Q %*% (lam * t(Q))
    if (residual) {
      total <- if (inherits(b, "bp_basis") && !is.null(b$total)) b$total else
        sum(lam)
      r <- max(total - sum(lam), 0) / max(1, p - length(j))
      Sh <- Sh + r * (diag(p) - Q %*% t(Q))
    }
    S <- S + weights[h] * Sh
  }
  S <- (S + t(S)) / 2
  e <- eigen(S, symmetric = TRUE)
  qq <- min(as.integer(q_out), p)
  out <- list(T = sign_fix(e$vectors[, seq_len(qq), drop = FALSE]),
              values = e$values[seq_len(qq)], mean = rep(0, p), q = qq,
              n = length(bases), p = p, method = "spectral")
  if (return_cov) list(T = out$T, S = S) else out
}

#' Pooled basis from a list of beat matrices
#'
#' Convenience wrapper that extracts the bases and pools them in one call.
#'
#' @param X_list A list of beat matrices.
#' @param q Number of components.
#' @param method Pooling method.
#' @param align Align before pooling.
#' @param iterations Number of GPA iterations.
#' @param weights Optional record weights.
#' @param ... Passed to [pool_bases()].
#' @return A \code{bp_pool} object.
#' @export
pool_from_beats <- function(X_list, q = 10L, method = "procrustes",
                            align = TRUE, iterations = 1L, weights = NULL,
                            ...) {
  bases <- extract_bases(X_list, q = q)
  pool_bases(bases, method = method, align = align,
             iterations = iterations, q = q, data = X_list,
             weights = weights, ...)
}

#' Leave-one-record-out pooled bases
#'
#' @param bases A list of bases.
#' @param data Optional list of beat matrices for column ordering.
#' @param method Pooling method.
#' @param align Align before pooling.
#' @param iterations Number of GPA iterations.
#' @param q Number of components.
#' @return A list with the full pooled basis \code{full} and a list
#'   \code{leave_one_out} of pooled bases obtained by omitting each record.
#' @export
leave_one_out_pool <- function(bases, data = NULL, method = "procrustes",
                               align = TRUE, iterations = 1L, q = NULL) {
  full <- pool_bases(bases, method = method, align = align,
                     iterations = iterations, q = q, data = data)
  H <- length(bases)
  loo <- vector("list", H)
  for (h in seq_len(H)) {
    b <- bases[-h]
    d <- if (!is.null(data)) data[-h] else NULL
    ref <- if (h == 1L) 1L else 1L
    loo[[h]] <- pool_bases(b, method = method, align = align,
                           iterations = iterations, q = q, data = d,
                           ref = ref)
  }
  list(full = full, leave_one_out = loo)
}

#' @export
as.matrix.bp_pool <- function(x, ...) x$T

#' @export
print.bp_pool <- function(x, ...) {
  cat("beatpool pooled basis\n")
  cat("  method:", x$method, "| aligned:", x$align,
      "| bases:", x$n_bases, "| components:", x$q, "\n")
  if (!is.null(x$values))
    cat("  explained variance (first 3):",
        paste(signif(utils::head(x$values, 3), 4), collapse = ", "), "\n")
  invisible(x)
}

#' @export
print.bp_basis <- function(x, ...) {
  cat("beatpool PCA basis\n")
  cat("  data:", x$n, "beats x", x$p, "samples | components:", x$q, "\n")
  cat("  explained variance (first 3):",
      paste(signif(utils::head(x$values, 3), 4), collapse = ", "), "\n")
  invisible(x)
}
