## Helpers shared by the article scripts. Each script is run from any
## directory; inputs and outputs live in BP_DATA and BP_OUT.

suppressPackageStartupMessages(library(beatpool))
data_dir <- Sys.getenv("BP_DATA", "data")
out_dir <- Sys.getenv("BP_OUT", "article-output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cores <- as.integer(Sys.getenv("BP_CORES", "4"))
if (.Platform$OS.type == "windows") cores <- 1L

CORPORA <- c(mitbih = "MIT-BIH", svdb = "SVDB", incart = "INCART", finance = "Finance")
out <- function(...) file.path(out_dir, ...)

## Centred beats of one record, their squared norms and the truncated
## energy-normalised moment M_h of Section 3.
record_stats <- function(X, q) {
  m <- energy_moment(X, q = q)
  Xc <- sweep(as.matrix(X), 2L, colMeans(X), "-")
  dn <- rowSums(Xc^2)
  list(Xc = Xc[dn > 0, , drop = FALSE], dn = dn[dn > 0], moment = m,
       M = m$T %*% (m$values * t(m$T)))
}
leading <- function(M, k)
  eigen((M + t(M)) / 2, symmetric = TRUE)$vectors[, seq_len(k), drop = FALSE]

## Mean per-beat relative error and ratio of traces of an orthonormal basis B.
beat_error <- function(Xc, dn, B) mean((dn - rowSums((Xc %*% B)^2)) / dn)
trace_error <- function(Xc, dn, B) 1 - sum((Xc %*% B)^2) / sum(dn)
choose_basis <- function(Bs, Xc, dn)
  which.min(vapply(Bs, function(B) beat_error(Xc, dn, B), numeric(1)))

## Balanced accuracy of two-cluster k-means and of cross-validated LDA.
kmeans_balanced <- function(cl, is_normal) {
  ab <- !is_normal
  tab <- table(cl, ab)
  if (nrow(tab) < 2L) return(NA_real_)
  cl_ab <- rownames(tab)[which.max(tab[, "TRUE"])]
  pred_ab <- cl == cl_ab
  (mean(pred_ab[ab]) + mean(!pred_ab[!ab])) / 2
}
lda_balanced <- function(Y, cls, folds = 5L, seed = 1L) {
  set.seed(seed)
  f <- integer(length(cls))
  for (g in unique(cls)) {
    i <- which(cls == g)
    f[i] <- sample(rep(seq_len(folds), length.out = length(i)))
  }
  acc <- numeric(folds)
  for (j in seq_len(folds)) {
    tr <- f != j
    fit <- tryCatch(MASS::lda(Y[tr, , drop = FALSE], cls[tr]), error = function(e) NULL)
    if (is.null(fit)) { acc[j] <- NA; next }
    pr <- predict(fit, Y[!tr, , drop = FALSE])$class
    tb <- table(pr, cls[!tr])
    if (nrow(tb) < 2L) { acc[j] <- NA; next }
    acc[j] <- mean(c(tb[2, 2] / sum(tb[, 2]), tb[1, 1] / sum(tb[, 1])))
  }
  mean(acc, na.rm = TRUE)
}
