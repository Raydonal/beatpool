#' All pooling operators at once
#'
#' Runs the operator families compared in the paper on a common set of local
#' bases and returns their \eqn{p \times q} operators, together with the number
#' of values each operator transmits per record. The function is the
#' programmatic counterpart of the comparison tables: the returned matrices can
#' be passed to [reconstruction_error()] or [subspace_distance()].
#'
#' @param X_list A list of beat matrices, one per record. Required unless
#'   \code{bases} is supplied.
#' @param bases Optional list of \code{bp_basis} objects. When \code{NULL},
#'   bases are extracted from \code{X_list}.
#' @param q Number of components.
#' @param weights Optional record weights.
#' @param residual Use residual imputation for the spectral summaries.
#' @return A list of class \code{bp_poolall}. Each element is a \eqn{p \times
#'   q} matrix with orthonormal columns, except the unaligned means, which are
#'   column-normalised but not orthogonal. The attribute \code{numbers} gives
#'   the communication cost of each operator per record.
#' @examples
#' sim <- simulate_beats(H = 6, n = 120, p = 20, q = 3, seed = 1)
#' all <- pool_all(sim$X, q = 3)
#' names(all)
#' attr(all, "numbers")
#' @export
pool_all <- function(X_list = NULL, bases = NULL, q = 10L, weights = NULL,
                     residual = FALSE) {
  if (is.null(bases)) {
    if (is.null(X_list)) {
      stop("supply X_list or bases", call. = FALSE)
    }
    if (!is.list(X_list)) X_list <- list(X_list)
    bases <- extract_bases(X_list, q = q)
  }
  if (inherits(bases, "bp_basis") || is.matrix(bases)) bases <- list(bases)
  if (!length(bases)) stop("no bases to pool", call. = FALSE)
  p <- nrow(.as_basis_matrix(bases[[1L]]))
  q <- min(as.integer(q), p, min(vapply(bases, function(b)
    ncol(.as_basis_matrix(b)), integer(1))))
  if (q < 1L) stop("q must be at least 1", call. = FALSE)

  out <- list(
    bar = pool_bases(bases, method = "bar", align = FALSE, q = q)$T,
    bara = pool_bases(bases, method = "bar", align = TRUE, q = q)$T,
    qr = pool_bases(bases, method = "qr", align = FALSE, q = q)$T,
    qra = pool_bases(bases, method = "qr", align = TRUE, q = q)$T,
    proc = pool_bases(bases, method = "procrustes", align = TRUE, q = q)$T,
    gpa = gpa_bases(bases, q = q)$T,
    hier = hierarchical_pool(bases, q = q)[[q]]$T,
    spec = pool_spectral(bases, q = q, q_out = q, weights = weights,
                         residual = residual)$T
  )
  if (q > 5L) {
    out$spec5 <- pool_spectral(bases, q = 5L, q_out = q, weights = weights,
                               residual = residual)$T
  }
  if (q > 3L) {
    out$spec3 <- pool_spectral(bases, q = 3L, q_out = q, weights = weights,
                               residual = residual)$T
  }
  numbers <- c(bar = p * q, bara = p * q, qr = p * q, qra = p * q,
               proc = p * q, gpa = p * q, hier = p * q,
               spec = q * (p + 1L))
  if (!is.null(out$spec5)) numbers["spec5"] <- 5L * (p + 1L)
  if (!is.null(out$spec3)) numbers["spec3"] <- 3L * (p + 1L)

  if (!is.null(X_list)) {
    if (!is.list(X_list)) X_list <- list(X_list)
    w <- if (is.null(weights)) vapply(X_list, nrow, integer(1)) else weights
    out$cov <- pool_covariance(X_list, q = q)$T
    out$covn <- pool_covariance(X_list, q = q, weights = w)$T
    numbers["cov"] <- p * (p + 1L) / 2
    numbers["covn"] <- p * (p + 1L) / 2
  }

  structure(out[order(names(out))], class = c("bp_poolall", "list"),
            numbers = numbers[order(names(out))], q = q, p = p)
}

#' @export
print.bp_poolall <- function(x, ...) {
  cat("beatpool operator collection\n")
  cat("  operators:", length(x), "| components:", attr(x, "q"),
      "| samples:", attr(x, "p"), "\n")
  num <- attr(x, "numbers")
  for (nm in names(x)) {
    cat(sprintf("  %-6s %5d values transmitted\n", nm, num[[nm]]))
  }
  invisible(x)
}

#' @export
summary.bp_pool <- function(object, ...) {
  data.frame(
    method = object$method,
    components = object$q,
    bases = object$n_bases,
    aligned = object$align,
    iterations = object$iterations_used,
    converged = object$converged,
    top_variance = if (is.null(object$values)) NA_real_ else object$values[1L],
    stringsAsFactors = FALSE
  )
}
