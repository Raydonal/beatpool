#' Simulate ECG-like heartbeat matrices
#'
#' Generates one or several records of beats from a Gaussian model with a
#' low-rank covariance. Heterogeneity across records is introduced by random
#' rotations of a common eigenbasis.
#'
#' @param H Number of records.
#' @param n Number of beats per record (a scalar or a vector of length
#'   \code{H}).
#' @param p Number of time samples.
#' @param q Number of signal components.
#' @param model \code{"homogeneous"}, \code{"heterogeneous"} or
#'   \code{"mixture"}.
#' @param epsilon Strength of the record-specific rotations.
#' @param mu Separation between the two class means in the mixture model.
#' @param pi Probability of the first class in the mixture model.
#' @param sigma Standard deviation of isotropic noise.
#' @param spectrum Optional vector of eigenvalues (length \code{q}).
#' @param basis Optional common \code{p x q} eigenbasis.
#' @param lambda_sd Standard deviation of log-normal record-specific
#'   eigenvalue perturbations.
#' @param seed Optional seed.
#' @return A list with the beat matrices \code{X}, the labels \code{labels}
#'   (\code{"normal"} and \code{"abnormal"} in the mixture model, \code{"normal"}
#'   for every beat in the other two, which have a single class),
#'   the normal indicators \code{normal}, the true record bases \code{Q_list},
#'   the common basis \code{Q0}, the eigenvalues \code{lambda}, and the
#'   simulation parameters.
#' @export
simulate_beats <- function(H = 20L, n = 200L, p = 61L, q = 10L,
                           model = c("homogeneous", "heterogeneous", "mixture"),
                           epsilon = 0, mu = 1, pi = 0.5, sigma = 0.2,
                           spectrum = NULL, basis = NULL, lambda_sd = 0,
                           seed = NULL) {
  model <- match.arg(model)
  if (!is.null(seed)) set.seed(seed)
  H <- as.integer(H); p <- as.integer(p); q <- as.integer(q)
  n <- rep(as.integer(n), length.out = H)
  if (is.null(spectrum)) spectrum <- exp(-seq(0, q - 1) / 2)
  if (length(spectrum) < q) spectrum <- rep(spectrum, length.out = q)
  lambda <- spectrum[seq_len(q)]
  if (is.null(basis)) basis <- dhmt_basis(p)[, seq_len(q), drop = FALSE]
  Q0 <- orthonormalize(basis, q)

  if (epsilon > 0 || model == "heterogeneous") {
    Q_list <- lapply(seq_len(H), function(h) {
      Omega <- matrix(stats::rnorm(p * p), p, p)
      Omega <- (Omega - t(Omega)) / sqrt(2 * p)
      R <- .expm(epsilon * Omega)
      orthonormalize(R %*% Q0, q)
    })
  } else {
    Q_list <- rep(list(Q0), H)
  }

  X <- vector("list", H)
  labels <- vector("list", H)
  normal <- vector("list", H)
  for (h in seq_len(H)) {
    lam <- lambda
    if (lambda_sd > 0) lam <- lam * exp(stats::rnorm(q, 0, lambda_sd))
    sd <- sqrt(lam)
    Z <- matrix(stats::rnorm(n[h] * q), n[h], q) %*% diag(sd, q)
    if (model == "mixture") {
      cl <- stats::rbinom(n[h], 1L, pi)
      dir <- Q0[, 1L]
      shift <- ifelse(cl == 1L, -mu, mu)
      Xh <- Z %*% t(Q_list[[h]]) + outer(shift, dir)
      labels[[h]] <- ifelse(cl == 1L, "abnormal", "normal")
      normal[[h]] <- cl == 0L
    } else {
      Xh <- Z %*% t(Q_list[[h]])
      labels[[h]] <- rep("normal", n[h])
      normal[[h]] <- rep(TRUE, n[h])
    }
    if (sigma > 0) {
      Xh <- Xh + matrix(stats::rnorm(n[h] * p, 0, sigma), n[h], p)
    }
    X[[h]] <- Xh
  }
  for (h in seq_len(H)) {
    attr(X[[h]], "record") <- paste0("sim", h)
    attr(X[[h]], "labels") <- labels[[h]]
    attr(X[[h]], "normal") <- normal[[h]]
  }
  list(X = X, labels = labels, normal = normal, Q_list = Q_list, Q0 = Q0,
       lambda = lambda, params = list(H = H, n = n, p = p, q = q,
                                      model = model, epsilon = epsilon,
                                      mu = mu, pi = pi, sigma = sigma,
                                      lambda_sd = lambda_sd))
}

#' Corrupt a basis with random signs and permutations
#'
#' Used to demonstrate the need for alignment before pooling.
#'
#' @param Q A basis matrix.
#' @param mode One of \code{"sign"}, \code{"permutation"} or \code{"both"}.
#' @param seed Optional seed.
#' @return A corrupted basis with the same span.
#' @export
corrupt_basis <- function(Q, mode = c("sign", "permutation", "both"),
                          seed = NULL) {
  mode <- match.arg(mode)
  if (!is.null(seed)) set.seed(seed)
  Q <- as.matrix(Q)
  out <- Q
  if (mode %in% c("sign", "both")) {
    s <- sample(c(-1, 1), ncol(Q), replace = TRUE)
    out <- sweep(out, 2L, s, "*")
  }
  if (mode %in% c("permutation", "both")) {
    out <- out[, sample.int(ncol(out)), drop = FALSE]
  }
  out
}

.expm <- function(A) {
  if (requireNamespace("expm", quietly = TRUE)) return(expm::expm(A))
  nrm <- norm(A, "O")
  s <- max(0L, ceiling(log2(max(1, nrm))))
  B <- A / 2^s
  E <- diag(nrow(A)); term <- E
  for (k in seq_len(20L)) {
    term <- term %*% B / k
    E <- E + term
  }
  for (i in seq_len(s)) E <- E %*% E
  E
}
