## Step 2. Held-out evaluation of every operator in the four corpora.
## One fold per record, or per patient in INCART. Writes, per corpus, the
## errors (heldout), the rate--distortion grid (rate), the basis chosen from
## all beats and from the first 50 (chosen), the normal-versus-ventricular
## task (downstream) and the grouping of all records (clusters).

source(file.path(Sys.getenv("BP_ARTICLE", system.file("article", package = "beatpool")), "helpers.R"))
suppressPackageStartupMessages(library(parallel))

QMAX <- 10L                       # eigenpairs transmitted per record, q
KS <- c(1L, 2L, 3L, 5L, 10L)      # decoding dimensions
KK <- c(3L, 10L)                  # dimensions with shared bases
NB <- 2:4                         # numbers of shared bases, C
N0 <- 50L                         # beats used to choose a basis early
STEPS <- 2^seq(0, -8, by = -0.5)  # quantiser steps
NREP <- 10L                       # balanced resamples in the task
NMAX <- 300L                      # beats per class in the task

run_corpus <- function(t) {
  X <- readRDS(out(sprintf("beats_%s.rds", t)))
  group <- if (t == "incart") vapply(X, attr, "", "patient") else names(X)
  task <- t != "finance"
  records <- names(X)
  H <- length(X)
  p <- ncol(X[[1]])
  nb <- vapply(X, nrow, integer(1))
  folds <- unique(group)
  bases <- lapply(X, function(x) pca_basis(x, q = QMAX))
  st <- lapply(X, record_stats, q = QMAX)
  DC <- dct_basis(p)

  one_fold <- function(fi) {
    held <- which(group == folds[fi])
    tr <- setdiff(seq_len(H), held)
    P <- list(
      spec = pool_spectral(bases[tr], q = QMAX)$T,
      covn = pool_covariance(X[tr], q = QMAX, weights = nb[tr])$T,
      proc = pool_bases(bases[tr], method = "procrustes", align = TRUE, ref = "medoid",
                        q = QMAX, data = X[tr])$T,
      enorm = leading(Reduce(`+`, lapply(st[tr], `[[`, "M")), QMAX),
      dct = DC)
    hier <- hierarchical_pool(bases[tr], q = QMAX, ref = "medoid", data = X[tr])
    shared <- list()
    for (k in KK) for (C in NB)
      shared[[paste(k, C)]] <- pool_shared(lapply(st[tr], `[[`, "moment"), k = k, C = C,
                                           seed = 1000L * fi + 10L * k + C)$T
    res <- list(); rd <- list(); ds <- list(); ch <- list()
    for (h in held) {
      z <- st[[h]]
      basis_of <- function(m, k) {
        if (m == "local") return(orthonormalize(bases[[h]]$T, k))
        if (m == "hier") return(orthonormalize(hier[[k]]$T, k))
        if (m == "dct") return(orthonormalize(DC, k))
        if (grepl("^multi", m)) {
          Bs <- shared[[paste(k, sub("multi", "", m))]]
          return(Bs[[choose_basis(Bs, z$Xc, z$dn)]])
        }
        orthonormalize(P[[m]], k)
      }
      add <- function(m, k, B)
        res[[length(res) + 1L]] <<- data.frame(
          corpus = CORPORA[[t]], record = records[h], method = m, k = k,
          beat = beat_error(z$Xc, z$dn, B), trace = trace_error(z$Xc, z$dn, B))
      for (k in KS) for (m in c("local", "spec", "covn", "proc", "hier", "enorm", "dct"))
        add(m, k, basis_of(m, k))
      n0 <- seq_len(min(N0, nrow(z$Xc)))
      for (k in KK) {
        for (C in NB) add(paste0("multi", C), k, basis_of(paste0("multi", C), k))
        Bs <- shared[[paste(k, 3L)]]
        i_all <- choose_basis(Bs, z$Xc, z$dn)
        i_n0 <- choose_basis(Bs, z$Xc[n0, , drop = FALSE], z$dn[n0])
        add("multi3_n0", k, Bs[[i_n0]])
        ch[[length(ch) + 1L]] <- data.frame(corpus = CORPORA[[t]], record = records[h],
                                            k = k, chosen = i_all, chosen_n0 = i_n0)
      }
      for (m in c("local", "spec", "enorm", "proc", "dct", "multi3")) {
        r <- rate_distortion(X[[h]], basis_of(m, 10L), k = 10L, steps = STEPS,
                             sample_bits = 16L, transmit_basis = FALSE)
        rd[[length(rd) + 1L]] <- data.frame(
          corpus = CORPORA[[t]], record = records[h], method = m, n = nrow(X[[h]]),
          nbasis = switch(m, dct = 0L, multi3 = 3L, 1L),
          step = r$step, rate0 = r$rate, D = r$distortion)
      }
      if (task) {
        lab <- attr(X[[h]], "labels")
        idx_n <- which(lab == "N")
        idx_v <- which(lab == "V")
        if (length(idx_n) >= 20L && length(idx_v) >= 20L) {
          nuse <- min(NMAX, length(idx_n), length(idx_v))
          meth <- c("local", "spec", "enorm", "proc", "multi3")
          Bs <- list()
          for (k in KK) for (m in meth) Bs[[paste(m, k)]] <- basis_of(m, k)
          for (b in seq_len(NREP)) {
            set.seed(50000L + 100L * h + b)
            sel <- c(sample(idx_n, nuse), sample(idx_v, nuse))
            Xs <- X[[h]][sel, , drop = FALSE]
            is_norm <- lab[sel] == "N"
            cls <- factor(ifelse(is_norm, "N", "V"), levels = c("N", "V"))
            for (k in KK) for (m in meth) {
              Y <- sweep(Xs, 2L, colMeans(Xs), "-") %*% Bs[[paste(m, k)]]
              cl <- stats::kmeans(Y, centers = 2L, nstart = 10L, iter.max = 100L)$cluster
              ds[[length(ds) + 1L]] <- data.frame(
                corpus = CORPORA[[t]], record = records[h], rep = b, method = m, k = k,
                bal = kmeans_balanced(cl, is_norm),
                lda_bal = lda_balanced(Y, cls, seed = 50000L + 100L * h + b))
            }
          }
        }
      }
    }
    list(res = do.call(rbind, res), rd = do.call(rbind, rd), ch = do.call(rbind, ch),
         ds = if (length(ds)) do.call(rbind, ds) else NULL)
  }

  t0 <- Sys.time()
  z <- mclapply(seq_along(folds), one_fold, mc.cores = cores, mc.preschedule = FALSE)
  bad <- vapply(z, inherits, logical(1), "try-error")
  if (any(bad)) stop(as.character(z[[which(bad)[1]]]))
  message(CORPORA[[t]], " done in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
  w <- function(key, file) {
    d <- do.call(rbind, lapply(z, `[[`, key))
    if (!is.null(d)) write.csv(d, out(sprintf("%s_%s.csv", t, file)), row.names = FALSE)
  }
  w("res", "heldout"); w("rd", "rate"); w("ch", "chosen"); w("ds", "downstream")

  full <- pool_shared(lapply(st, `[[`, "moment"), k = 3L, C = 3L, seed = 7L)
  write.csv(data.frame(
    record = records, group = group, cluster = full$group, beats = nb,
    normal = vapply(X, function(x) mean(attr(x, "normal")), numeric(1))),
    out(sprintf("%s_clusters.csv", t)), row.names = FALSE)
}

for (t in names(CORPORA)) run_corpus(t)
