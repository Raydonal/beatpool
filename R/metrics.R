#' Relative reconstruction error
#'
#' Projects each beat onto the subspace spanned by the first \code{k} columns
#' of a basis and returns the relative squared error.
#'
#' @param X A beat matrix.
#' @param T A basis matrix.
#' @param k Number of leading components.
#' @param center Center each column of \code{X} before projecting.
#' @param per_beat Return one value per beat instead of the average.
#' @return A numeric scalar or vector.
#' @export
reconstruction_error <- function(X, T, k = ncol(T), center = TRUE,
                                 per_beat = FALSE) {
  X <- as.matrix(X)
  Q <- orthonormalize(T, k = k)
  Xc <- if (center) sweep(X, 2L, colMeans(X), "-") else X
  P <- Xc %*% Q %*% t(Q)
  num <- rowSums((Xc - P)^2)
  den <- rowSums(Xc^2)
  err <- ifelse(den > 0, num / den, 0)
  if (per_beat) err else mean(err)
}

#' Reconstruction error curve
#'
#' @param X A beat matrix.
#' @param T A basis matrix.
#' @param k Integer vector of numbers of components.
#' @param center Center before projecting.
#' @return A data frame with columns \code{k} and \code{error}.
#' @export
reconstruction_curve <- function(X, T, k = seq_len(ncol(T)), center = TRUE) {
  data.frame(k = k,
             error = vapply(k, function(kk)
               reconstruction_error(X, T, kk, center = center), numeric(1)))
}

#' Compress and reconstruct beats
#'
#' @param X A beat matrix.
#' @param T A basis matrix.
#' @param k Number of components.
#' @param center Center before projecting.
#' @return A list with the coefficients \code{Y}, the reconstructed beats
#'   \code{X_hat}, the per-beat relative errors \code{error} and the mean
#'   relative error \code{mean_error}.
#' @export
compress_beats <- function(X, T, k = ncol(T), center = TRUE) {
  X <- as.matrix(X)
  Q <- orthonormalize(T, k = k)
  mu <- if (center) colMeans(X) else rep(0, ncol(X))
  Xc <- sweep(X, 2L, mu, "-")
  Y <- Xc %*% Q
  Xh <- Y %*% t(Q)
  err <- rowSums((Xc - Xh)^2) / pmax(rowSums(Xc^2), .Machine$double.eps)
  list(Y = Y, X_hat = sweep(Xh, 2L, mu, "+"), error = err,
       mean_error = mean(err), k = k, T = Q, mean = mu)
}

#' Rate--distortion curve of a quantised spectral codec
#'
#' Quantises the coefficients of a linear projection with a uniform scalar
#' quantiser whose step is a multiple of the coefficient standard deviation,
#' entropy-codes the integer levels, and reports the mean relative distortion
#' and the rate in bits per beat. The basis and the mean are charged once per
#' data holder and amortised over \code{n_amortize} beats, which is the cost
#' model used in the paper.
#'
#' @param X A beat matrix, one beat per row.
#' @param T A basis matrix.
#' @param k Number of components.
#' @param steps Multipliers of the coefficient standard deviation; a larger
#'   multiplier gives a coarser quantiser, a lower rate and a larger
#'   distortion.
#' @param center Center before projecting.
#' @param sample_bits Bits per raw sample in the reference codec, used only
#'   for the compression ratio.
#' @param basis_bits Bits per entry of the transmitted basis.
#' @param mean_bits Bits per entry of the transmitted mean.
#' @param n_amortize Number of beats over which the basis and the mean are
#'   amortised; defaults to the number of beats in \code{X}.
#' @param transmit_basis Whether the basis is transmitted. Fixed transforms
#'   such as the DCT are known to the decoder and cost nothing.
#' @return A data frame with the step multiplier \code{step}, the rate in bits
#'   per beat, the mean relative distortion in percent and the compression
#'   ratio against the raw reference.
#' @export
rate_distortion <- function(X, T, k = ncol(T),
                            steps = 2^seq(0, -8, by = -0.5),
                            center = TRUE, sample_bits = 11L,
                            basis_bits = 16L, mean_bits = 16L,
                            n_amortize = nrow(X), transmit_basis = TRUE) {
  X <- as.matrix(X)
  p <- ncol(X)
  Q <- orthonormalize(T, k = k)
  mu <- if (center) colMeans(X) else rep(0, p)
  Xc <- sweep(X, 2L, mu, "-")
  Y <- Xc %*% Q
  s <- apply(Y, 2L, stats::sd)
  s[!is.finite(s) | s <= 0] <- 1
  ref <- sample_bits * p
  overhead <- p * mean_bits + ncol(Q) * 16 +
    if (transmit_basis) ncol(Q) * p * basis_bits else 0
  out <- lapply(steps, function(a) {
    step <- a * s
    idx <- round(sweep(Y, 2L, step, "/"))
    Yh <- sweep(idx, 2L, step, "*")
    err <- rowSums((Xc - Yh %*% t(Q))^2) /
      pmax(rowSums(Xc^2), .Machine$double.eps)
    coef_bits <- sum(vapply(seq_len(ncol(idx)), function(j)
      nrow(idx) * .entropy_bits(idx[, j]), numeric(1)))
    rate <- (coef_bits + overhead) / n_amortize
    data.frame(step = a, rate = rate, distortion = 100 * mean(err),
               compression_ratio = ref / rate)
  })
  do.call(rbind, out)
}

.entropy_bits <- function(idx) {
  pr <- table(idx) / length(idx)
  -sum(pr * log2(pr))
}

#' Adjusted Rand index
#'
#' @param x Cluster labels.
#' @param y True labels.
#' @return The adjusted Rand index, or \code{NA} when it is undefined.
#' @export
ari <- function(x, y) {
  x <- as.character(x); y <- as.character(y)
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 2L) return(NA_real_)
  tab <- table(x, y)
  comb2 <- function(z) sum(z * (z - 1) / 2)
  idx <- comb2(tab)
  a <- sum(comb2(rowSums(tab)))
  b <- sum(comb2(colSums(tab)))
  expected <- a * b / (n * (n - 1) / 2)
  maxidx <- (a + b) / 2
  den <- maxidx - expected
  if (!is.finite(den) || abs(den) < .Machine$double.eps) return(NA_real_)
  (idx - expected) / den
}

#' Local cluster consistency
#'
#' For each beat, the proportion of correct assignments among the other beats
#' in the same cluster must exceed one half. This is the criterion used in the
#' original ECG analysis.
#'
#' @param cluster Cluster labels.
#' @param truth True labels.
#' @return The mean proportion of locally consistent beats.
#' @export
prop_correct <- function(cluster, truth) {
  cluster <- as.character(cluster); truth <- as.character(truth)
  n <- length(cluster)
  if (length(truth) != n) stop("vectors must have equal length", call. = FALSE)
  correct <- logical(n)
  for (i in seq_len(n)) {
    same <- which(cluster == cluster[i])
    same <- same[same != i]
    if (!length(same)) { correct[i] <- FALSE; next }
    correct[i] <- mean(truth[i] == truth[same]) > 0.5
  }
  mean(correct)
}

#' Cluster purity
#'
#' @param cluster Cluster labels.
#' @param truth True labels.
#' @return The weighted average of the largest class proportion in each
#'   cluster.
#' @export
purity <- function(cluster, truth) {
  tab <- table(cluster, truth)
  sum(apply(tab, 1L, max)) / sum(tab)
}

#' Cluster quality metrics
#'
#' @param cluster Cluster labels.
#' @param truth True labels.
#' @return A one-row data frame with \code{prop_correct}, \code{ari} and
#'   \code{purity}.
#' @export
cluster_metrics <- function(cluster, truth) {
  data.frame(prop_correct = prop_correct(cluster, truth),
             ari = ari(cluster, truth),
             purity = purity(cluster, truth))
}

#' Cluster beats in a projected basis
#'
#' @param X A beat matrix.
#' @param T A basis matrix.
#' @param k Number of leading components used as features.
#' @param K Number of clusters.
#' @param truth Optional true labels used to compute quality metrics.
#' @param center Center before projecting.
#' @param nstart Number of random starts.
#' @param iter.max Maximum number of k-means iterations.
#' @param seed Optional seed.
#' @return A list with the cluster labels, the fitted \code{kmeans} object, the
#'   features \code{Y} and, when \code{truth} is supplied, \code{metrics}.
#' @export
cluster_beats <- function(X, T, k = 3L, K = 2L, truth = NULL, center = TRUE,
                          nstart = 10L, iter.max = 100L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Q <- orthonormalize(T, k = k)
  Xc <- if (center) sweep(X, 2L, colMeans(X), "-") else as.matrix(X)
  Y <- Xc %*% Q
  fit <- stats::kmeans(Y, centers = K, nstart = nstart, iter.max = iter.max)
  out <- list(cluster = fit$cluster, fit = fit, Y = Y, k = k, K = K)
  if (!is.null(truth)) out$metrics <- cluster_metrics(fit$cluster, truth)
  out
}

#' Jackknife stability of a pooled basis
#'
#' @param loo A list returned by [leave_one_out_pool()].
#' @param type Distance used between the full and leave-one-out bases.
#' @return A data frame with the relative Frobenius distance and the subspace
#'   distance for each omitted record.
#' @export
jackknife_stability <- function(loo, type = "rms") {
  full <- loo$full$T
  data.frame(
    record = seq_along(loo$leave_one_out),
    relative = vapply(loo$leave_one_out, function(z)
      norm(z$T - full, "F") / norm(full, "F"), numeric(1)),
    subspace = vapply(loo$leave_one_out, function(z)
      subspace_distance(z$T, full, type = type), numeric(1))
  )
}

#' Percentile bootstrap confidence interval
#'
#' @param x Numeric vector.
#' @param B Number of bootstrap replicates.
#' @param conf Confidence level.
#' @param statistic Function returning the statistic of interest.
#' @param seed Optional seed.
#' @return A named numeric vector with the lower limit, estimate and upper
#'   limit.
#' @export
bootstrap_ci <- function(x, B = 2000L, conf = 0.95, statistic = mean,
                         seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x <- x[is.finite(x)]
  n <- length(x)
  if (!n) return(c(lower = NA, estimate = NA, upper = NA))
  est <- statistic(x)
  bs <- vapply(seq_len(B), function(b)
    statistic(x[sample.int(n, n, replace = TRUE)]), numeric(1))
  alpha <- (1 - conf) / 2
  c(lower = unname(stats::quantile(bs, alpha)),
    estimate = est,
    upper = unname(stats::quantile(bs, 1 - alpha)))
}

#' Paired bootstrap comparison
#'
#' @param x Numeric vector for method one.
#' @param y Numeric vector for method two, paired with \code{x}.
#' @param B Number of bootstrap replicates.
#' @param conf Confidence level.
#' @param seed Optional seed.
#' @return A list with the observed mean difference, its confidence interval
#'   and a two-sided bootstrap p-value.
#' @export
paired_bootstrap <- function(x, y, B = 2000L, conf = 0.95, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  d <- x - y
  ok <- is.finite(d)
  d <- d[ok]
  n <- length(d)
  if (!n) return(list(estimate = NA, lower = NA, upper = NA, p_value = NA))
  bs <- vapply(seq_len(B), function(b) mean(d[sample.int(n, n, TRUE)]),
               numeric(1))
  alpha <- (1 - conf) / 2
  list(estimate = mean(d),
       lower = unname(stats::quantile(bs, alpha)),
       upper = unname(stats::quantile(bs, 1 - alpha)),
       p_value = min(1, 2 * min(mean(bs <= 0), mean(bs >= 0))))
}

#' Compare several methods across records
#'
#' @param values A data frame with one row per record and one column per
#'   method.
#' @param reference Name of the reference method.
#' @return A data frame with the mean, bootstrap interval and paired
#'   comparison against the reference for each method.
#' @export
compare_methods <- function(values, reference = colnames(values)[1]) {
  methods <- colnames(values)
  ref <- values[[reference]]
  do.call(rbind, lapply(methods, function(m) {
    x <- values[[m]]
    ci <- bootstrap_ci(x)
    if (m == reference) {
      data.frame(method = m, mean = ci["estimate"], lower = ci["lower"],
                 upper = ci["upper"], diff = 0, diff_lower = 0,
                 diff_upper = 0, p_value = NA_real_)
    } else {
      pb <- paired_bootstrap(x, ref)
      data.frame(method = m, mean = ci["estimate"], lower = ci["lower"],
                 upper = ci["upper"], diff = pb$estimate,
                 diff_lower = pb$lower, diff_upper = pb$upper,
                 p_value = pb$p_value)
    }
  }))
}
