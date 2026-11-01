#' Energy-normalised moment of one record
#'
#' Centres the beats, divides each one by its norm and returns the leading
#' eigenpairs of the resulting moment
#' \eqn{R_h = n_h^{-1} \sum_i x_{hi} x_{hi}^T / \|x_{hi}\|^2}. For a matrix
#' \eqn{T} with orthonormal columns, the mean per-beat relative error of the
#' record is \eqn{1 - tr(T^T R_h T)}, so that the moment targets that error
#' directly. The sample covariance targets the ratio of traces instead.
#'
#' @param X A beat matrix with beats in rows.
#' @param q Number of eigenpairs retained.
#' @return An object of class \code{bp_moment} with the eigenvectors \code{T},
#'   the eigenvalues \code{values}, the unit-trace tail \code{tail} and the
#'   number of beats \code{n}.
#' @examples
#' sim <- simulate_beats(H = 2, n = 50, p = 20, q = 3, seed = 1)
#' m <- energy_moment(sim$X[[1]], q = 5)
#' m$tail
#' @export
energy_moment <- function(X, q = 10L) {
  X <- as.matrix(X)
  Xc <- sweep(X, 2L, colMeans(X), "-")
  dn <- rowSums(Xc^2)
  Xc <- Xc[dn > 0, , drop = FALSE]
  dn <- dn[dn > 0]
  R <- crossprod(Xc / sqrt(dn)) / nrow(Xc)
  e <- eigen((R + t(R)) / 2, symmetric = TRUE)
  q <- min(as.integer(q), ncol(X))
  j <- seq_len(q)
  lam <- pmax(e$values, 0)
  structure(list(T = e$vectors[, j, drop = FALSE], values = lam[j],
                 tail = sum(lam[-j]), n = nrow(Xc), p = ncol(X), q = q),
            class = "bp_moment")
}

.moment_matrix <- function(m) m$T %*% (m$values * t(m$T))

.leading <- function(M, k) {
  eigen((M + t(M)) / 2, symmetric = TRUE)$vectors[, seq_len(k), drop = FALSE]
}

#' Pool energy-normalised moments
#'
#' Diagonalises the average of the rank-q moments transmitted by the records.
#' The cost per record is the same as that of [pool_spectral()].
#'
#' @param moments A list of \code{bp_moment} objects.
#' @param q_out Number of components returned.
#' @param weights Optional record weights.
#' @return A \code{bp_basis} object.
#' @examples
#' sim <- simulate_beats(H = 4, n = 50, p = 20, q = 3, seed = 2)
#' m <- lapply(sim$X, energy_moment, q = 5)
#' pool_energy(m, q_out = 3)$T[1:2, ]
#' @export
pool_energy <- function(moments, q_out = moments[[1]]$q, weights = NULL) {
  if (is.null(weights)) weights <- rep(1, length(moments))
  weights <- weights / sum(weights)
  S <- Reduce(`+`, Map(function(m, w) w * .moment_matrix(m), moments, weights))
  e <- eigen((S + t(S)) / 2, symmetric = TRUE)
  j <- seq_len(as.integer(q_out))
  list(T = sign_fix(e$vectors[, j, drop = FALSE]), values = e$values[j],
       mean = rep(0, nrow(S)), q = length(j), n = length(moments), p = nrow(S),
       method = "energy")
}

#' Several shared bases from transmitted moments
#'
#' Partitions the records into \code{C} groups and returns one basis per group.
#' The receiver alternates between the leading \code{k} eigenvectors of each
#' group sum and the assignment that retains the most of each moment, from
#' \code{starts} random assignments. It uses only what the records have already
#' transmitted, so no further communication from them is needed.
#'
#' @param moments A list of \code{bp_moment} objects.
#' @param k Dimension of each basis.
#' @param C Number of shared bases.
#' @param starts Number of random initial assignments.
#' @param iters Maximum number of alternations per start.
#' @param seed Optional seed. The random state of the session is restored.
#' @return An object of class \code{bp_shared} with the list of bases \code{T},
#'   the assignment \code{group} and the objective \code{value}.
#' @examples
#' sim <- simulate_beats(H = 6, n = 50, p = 20, q = 3, seed = 3)
#' m <- lapply(sim$X, energy_moment, q = 5)
#' fit <- pool_shared(m, k = 2, C = 2, starts = 2, seed = 1)
#' fit$group
#' @export
pool_shared <- function(moments, k, C = 3L, starts = 10L, iters = 50L,
                        seed = NULL) {
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir = globalenv()))
      get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()))
    set.seed(seed)
  }
  Ms <- lapply(moments, .moment_matrix)
  m <- length(Ms)
  if (C > m) stop("`C` cannot exceed the number of records.", call. = FALSE)
  best <- NULL
  for (s in seq_len(starts)) {
    g <- sample(rep_len(seq_len(C), m))
    for (it in seq_len(iters)) {
      Ts <- lapply(seq_len(C), function(c) {
        if (!any(g == c)) g[sample.int(m, 1L)] <- c
        .leading(Reduce(`+`, Ms[g == c]), k)
      })
      sc <- vapply(Ts, function(Tc)
        vapply(Ms, function(M) sum(Tc * (M %*% Tc)), numeric(1)), numeric(m))
      ng <- max.col(matrix(sc, m, C), ties.method = "first")
      if (all(ng == g)) break
      g <- ng
    }
    v <- sum(matrix(sc, m, C)[cbind(seq_len(m), g)])
    if (is.null(best) || v > best$value)
      best <- list(value = v, T = Ts, group = g)
  }
  structure(c(best, list(k = k, C = C)), class = "bp_shared")
}

#' Choose a shared basis for a record
#'
#' Returns the index of the basis with the lowest mean per-beat relative error
#' on the given beats. A compressor transmits this index, of \code{log2(C)} bits.
#'
#' @param fit A \code{bp_shared} object.
#' @param X A beat matrix.
#' @return An integer index, with the errors in the attribute \code{errors}.
#' @examples
#' sim <- simulate_beats(H = 6, n = 50, p = 20, q = 3, seed = 3)
#' m <- lapply(sim$X, energy_moment, q = 5)
#' fit <- pool_shared(m, k = 2, C = 2, starts = 2, seed = 1)
#' select_shared(fit, sim$X[[1]])
#' @export
select_shared <- function(fit, X) {
  e <- vapply(fit$T, function(Q) reconstruction_error(X, Q, k = ncol(Q)),
              numeric(1))
  structure(which.min(e), errors = e)
}

#' @export
print.bp_shared <- function(x, ...) {
  cat("<bp_shared> ", x$C, " bases of dimension ", x$k, ", group sizes ",
      paste(tabulate(x$group, x$C), collapse = "/"), "\n", sep = "")
  invisible(x)
}
