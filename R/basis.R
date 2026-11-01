#' Fix the sign of the columns of a basis
#'
#' Multiplies each column by -1 when its largest absolute entry is negative.
#' This makes the basis reproducible without changing its span.
#'
#' @param T A matrix with orthonormal (or arbitrary) columns.
#' @return The sign-adjusted matrix.
#' @export
sign_fix <- function(T) {
  T <- as.matrix(T)
  if (!ncol(T)) return(T)
  for (j in seq_len(ncol(T))) {
    v <- T[, j]
    k <- which.max(abs(v))
    if (length(k) && v[k] < 0) T[, j] <- -v
  }
  T
}

.check_numeric_matrix <- function(X, what = "X", min_rows = 1L) {
  if (!is.matrix(X)) {
    if (is.data.frame(X)) X <- as.matrix(X) else
      stop(what, " must be a numeric matrix", call. = FALSE)
  }
  if (!is.numeric(X)) stop(what, " must be numeric", call. = FALSE)
  if (nrow(X) < min_rows) {
    stop(what, " must have at least ", min_rows, " rows", call. = FALSE)
  }
  if (ncol(X) < 1L) stop(what, " must have at least one column", call. = FALSE)
  if (any(!is.finite(X))) stop(what, " contains non-finite values", call. = FALSE)
  X
}

#' Principal component basis of a beat matrix
#'
#' @param X A numeric matrix with beats in rows and time samples in columns.
#' @param q Number of components to keep.
#' @param center Center each time sample before the decomposition.
#' @param fix_sign Apply [sign_fix()] to the returned basis.
#' @param scores Return the projected scores as well.
#' @return A list with \code{T} (the \code{p x q} basis), \code{values}
#'   (explained variances), \code{mean} (column means), \code{total} (total
#'   variance), \code{q} and, when requested, \code{scores}.
#' @export
pca_basis <- function(X, q = 10L, center = TRUE, fix_sign = TRUE,
                      scores = FALSE) {
  X <- .check_numeric_matrix(X, "X", min_rows = 1L)
  p <- ncol(X)
  q <- min(as.integer(q), p, nrow(X))
  mu <- if (center) colMeans(X) else rep(0, p)
  Xc <- sweep(X, 2L, mu, "-")
  s <- svd(Xc, nu = 0L, nv = q)
  T <- s$v[, seq_len(q), drop = FALSE]
  if (fix_sign) T <- sign_fix(T)
  values <- (s$d[seq_len(q)]^2) / max(1, nrow(X) - 1L)
  out <- list(T = T, values = values, mean = mu, q = q,
              n = nrow(X), p = p, method = "pca",
              total = sum(Xc^2) / max(1, nrow(X) - 1L))
  if (scores) out$scores <- Xc %*% T
  class(out) <- c("bp_basis", "list")
  out
}

#' Discrete Hermite transform basis
#'
#' Builds the eigenvectors of the tridiagonal matrix used by Mugler's discrete
#' Hermite transform. The construction follows the one used in the original
#' ECG analysis.
#'
#' @param p Number of time samples.
#' @param b Bandwidth parameter.
#' @return A \code{p x p} orthogonal matrix.
#' @export
dhmt_basis <- function(p, b = 1) {
  p <- as.integer(p)
  Tb <- matrix(0, p, p)
  tau <- 1 / (p * b^2)
  for (mu in 0:(p - 1L)) {
    k <- mu + 1L
    Tb[k, k] <- (-2) * cos(pi * p * tau) *
      sin(pi * mu * tau) * sin(pi * (p - mu - 1L) * tau)
    if (mu > 0L) {
      aux <- sin(pi * mu * tau) * sin(pi * (p - mu) * tau)
      Tb[k, k - 1L] <- aux
      Tb[k - 1L, k] <- aux
    }
  }
  sign_fix(eigen(Tb, symmetric = TRUE)$vectors)
}

#' Discrete cosine transform basis (DCT-II)
#'
#' @param p Number of time samples.
#' @return A \code{p x p} orthogonal matrix.
#' @export
dct_basis <- function(p) {
  p <- as.integer(p)
  j <- 0:(p - 1L)
  k <- 0:(p - 1L)
  T <- outer(j, k, function(j, k) cos(pi * (j + 0.5) * k / p)) *
    sqrt(2 / p)
  T[, 1L] <- T[, 1L] / sqrt(2)
  sign_fix(T)
}

#' Random orthogonal basis
#'
#' @param p Number of rows.
#' @param q Number of columns.
#' @param seed Optional seed.
#' @return A \code{p x q} matrix with orthonormal columns.
#' @export
random_basis <- function(p, q = p, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  q <- min(as.integer(q), as.integer(p))
  qr.Q(qr(matrix(stats::rnorm(p * q), p, q)))
}

#' Average explained variance of each basis column
#'
#' @param T A basis matrix.
#' @param X_list A list of beat matrices used to score the columns.
#' @param weights Optional record weights.
#' @return A numeric vector with one value per column of \code{T}.
#' @export
basis_variance <- function(T, X_list, weights = NULL) {
  T <- as.matrix(T)
  if (!is.list(X_list)) X_list <- list(X_list)
  if (is.null(weights)) weights <- rep(1, length(X_list))
  weights <- weights / sum(weights)
  v <- numeric(ncol(T))
  for (h in seq_along(X_list)) {
    X <- as.matrix(X_list[[h]])
    Xc <- sweep(X, 2L, colMeans(X), "-")
    Z <- Xc %*% T
    v <- v + weights[h] * colSums(Z^2) / max(1, nrow(Z) - 1L)
  }
  v
}

#' Order the columns of a basis by explained variance
#'
#' @param T A basis matrix.
#' @param X_list A list of beat matrices.
#' @param weights Optional record weights.
#' @param fix_sign Apply [sign_fix()] after reordering.
#' @return A list with the reordered basis \code{T} and the variances
#'   \code{values}.
#' @export
order_basis <- function(T, X_list, weights = NULL, fix_sign = TRUE) {
  v <- basis_variance(T, X_list, weights)
  o <- order(v, decreasing = TRUE)
  T <- as.matrix(T)[, o, drop = FALSE]
  if (fix_sign) T <- sign_fix(T)
  list(T = T, values = v[o], order = o)
}

#' Orthonormalise the first \code{k} columns of a basis
#'
#' Returns an orthonormal basis for the span of the selected columns. This
#' puts non-orthogonal pooling rules such as the column-normalised mean on the
#' same footing as orthogonal ones when computing projections or clustering
#' features.
#'
#' @param T A basis matrix.
#' @param k Number of leading columns to use.
#' @param fix_sign Apply [sign_fix()] to the result.
#' @return A \code{p x k} matrix with orthonormal columns.
#' @export
orthonormalize <- function(T, k = ncol(T), fix_sign = TRUE) {
  T <- as.matrix(T)[, seq_len(min(k, ncol(T))), drop = FALSE]
  if (ncol(T) == 1L) {
    q <- T / sqrt(sum(T^2))
  } else {
    q <- qr.Q(qr(T))
  }
  if (fix_sign) q <- sign_fix(q)
  q
}

#' Align a basis to a reference basis
#'
#' Finds the orthogonal matrix \code{O} that maximises
#' \code{trace(t(ref) \%*\% Q \%*\% O)}. This removes the sign, permutation and
#' rotation indeterminacy of principal component bases before they are pooled.
#'
#' @param Q A basis to align.
#' @param ref The reference basis.
#' @return A list with the aligned basis \code{Q}, the rotation \code{O} and
#'   the achieved correlation \code{correlation}.
#' @export
align_basis <- function(Q, ref) {
  Q <- as.matrix(Q); ref <- as.matrix(ref)
  q <- min(ncol(Q), ncol(ref))
  Q <- Q[, seq_len(q), drop = FALSE]
  ref <- ref[, seq_len(q), drop = FALSE]
  M <- crossprod(ref, Q)
  s <- svd(M)
  O <- s$v %*% t(s$u)
  list(Q = Q %*% O, O = O, correlation = sum(s$d))
}

#' Principal angles between two subspaces
#'
#' @param T1 A matrix whose columns span the first subspace.
#' @param T2 A matrix whose columns span the second subspace.
#' @return A numeric vector of principal angles in radians, in decreasing order.
#' @export
principal_angles <- function(T1, T2) {
  T1 <- orthonormalize(T1); T2 <- orthonormalize(T2)
  q <- min(ncol(T1), ncol(T2))
  s <- svd(crossprod(T1[, seq_len(q), drop = FALSE],
                     T2[, seq_len(q), drop = FALSE]))$d
  acos(pmin(1, pmax(0, s)))
}

#' Distance between two subspaces
#'
#' @param T1 A matrix whose columns span the first subspace.
#' @param T2 A matrix whose columns span the second subspace.
#' @param type \code{"rms"} for the root mean squared sine of the principal
#'   angles, \code{"max"} for the largest angle, or \code{"frobenius"} for
#'   \code{||T1 T1' - T2 T2'||_F / sqrt(2 q)}.
#' @return A non-negative scalar.
#' @export
subspace_distance <- function(T1, T2, type = c("rms", "max", "frobenius")) {
  type <- match.arg(type)
  T1 <- orthonormalize(T1); T2 <- orthonormalize(T2)
  q <- min(ncol(T1), ncol(T2))
  T1 <- T1[, seq_len(q), drop = FALSE]
  T2 <- T2[, seq_len(q), drop = FALSE]
  if (type == "frobenius") {
    return(norm(T1 %*% t(T1) - T2 %*% t(T2), "F") / sqrt(2 * q))
  }
  s <- svd(crossprod(T1, T2))$d
  if (type == "max") {
    acos(pmin(1, pmax(0, s[1L])))
  } else {
    sqrt(mean(pmax(0, 1 - s^2)))
  }
}

#' Convert a basis to a tidy data frame
#'
#' @param T A basis matrix.
#' @return A data frame with columns \code{time}, \code{component} and
#'   \code{value}.
#' @export
basis_frame <- function(T) {
  T <- as.matrix(T)
  data.frame(
    time = rep(seq_len(nrow(T)), times = ncol(T)),
    component = rep(seq_len(ncol(T)), each = nrow(T)),
    value = as.vector(T)
  )
}

#' Complete a partial basis to a full orthogonal matrix
#'
#' @param T A matrix with orthonormal columns.
#' @return A square orthogonal matrix whose first columns are \code{T}.
#' @export
complete_basis <- function(T) {
  T <- orthonormalize(T, ncol(T))
  qr.Q(qr(T), complete = TRUE)
}
