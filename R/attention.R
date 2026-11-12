## Attention over transmitted moments.
##
## Keys are the projectors onto the leading k eigenvectors of each transmitted
## moment M_h, values are the M_h themselves, and the query is the energy
## normalised moment R_0 of beats of a new record. The internal engine below is
## shared by the exported functions; it works on a prepared set of keys and
## values so that leave-one-out searches only subset columns.

.att_prep <- function(moments, ks) {
  Ms <- lapply(moments, .moment_matrix)
  p <- nrow(Ms[[1]])
  V <- vapply(Ms, as.numeric, numeric(p * p))
  ks <- sort(unique(as.integer(ks)))
  keys <- lapply(ks, function(k) lapply(moments, function(m)
    m$T[, seq_len(min(k, m$q)), drop = FALSE]))
  names(keys) <- ks
  .att_pop(list(V = V, p = p, H = length(Ms), keys = keys, ks = ks,
                vnorm = sqrt(colSums(V^2))))
}

## Mean of the values, scaled to unit trace, and its leading eigenvectors.
.att_pop <- function(prep) {
  M <- matrix(rowMeans(prep$V), prep$p, prep$p)
  prep$Mbar <- M / sum(diag(M))
  e <- eigen((prep$Mbar + t(prep$Mbar)) / 2, symmetric = TRUE)$vectors
  prep$Mlead <- lapply(prep$ks, function(k) e[, seq_len(k), drop = FALSE])
  names(prep$Mlead) <- prep$ks
  prep
}

.att_subset <- function(prep, idx) {
  prep$V <- prep$V[, idx, drop = FALSE]
  prep$vnorm <- prep$vnorm[idx]
  prep$keys <- lapply(prep$keys, `[`, idx)
  prep$H <- length(idx)
  .att_pop(prep)
}

.key_scores <- function(R0, keys, k, score, P0 = NULL) {
  switch(score,
    energy = vapply(keys, function(Tk) sum(Tk * (R0 %*% Tk)), numeric(1)),
    projection = vapply(keys, function(Tk) sum(crossprod(P0, Tk)^2) / k,
                        numeric(1)),
    stop("Unknown score.", call. = FALSE))
}

.att_scores <- function(R0, prep, k, score) {
  if (score == "cosine") {
    r <- as.numeric(R0)
    return(as.numeric(crossprod(r, prep$V)) / (sqrt(sum(r^2)) * prep$vnorm))
  }
  P0 <- if (score == "projection") .leading(R0, k) else NULL
  .key_scores(R0, prep$keys[[as.character(k)]], k, score, P0)
}

.sink_score <- function(R0, prep, k, score) {
  if (score == "cosine")
    return(sum(R0 * prep$Mbar) / (sqrt(sum(R0^2)) * sqrt(sum(prep$Mbar^2))))
  P0 <- if (score == "projection") .leading(R0, k) else NULL
  .key_scores(R0, list(prep$Mlead[[as.character(k)]]), k, score, P0)
}

.sparsemax <- function(z) {
  s <- sort(z, decreasing = TRUE)
  cs <- cumsum(s)
  j <- seq_along(s)
  kk <- max(j[1 + j * s > cs])
  pmax(z - (cs[kk] - 1) / kk, 0)
}

.softmax <- function(s, tau) {
  if (!is.finite(tau)) return(rep(1 / length(s), length(s)))
  if (tau <= 0) {
    w <- as.numeric(s >= max(s) - 1e-12)
    return(w / sum(w))
  }
  z <- exp((s - max(s)) / tau)
  z / sum(z)
}

.normalise <- function(s, tau, normalize) {
  if (normalize == "softmax" || !is.finite(tau) || tau <= 0)
    return(.softmax(s, tau))
  if (normalize == "sparsemax") return(.sparsemax(s / tau))
  stop("Unknown normalisation.", call. = FALSE)
}

## Attended moment for one query. Each head is scaled to unit trace before
## heads are averaged and before the contrast and self terms are applied.
.att_core <- function(R0, prep, k, tau, score = "energy", normalize = "softmax",
                      top = NULL, sink = NULL, contrast = 0, self = 0,
                      heads = NULL) {
  hk <- if (is.null(heads)) as.integer(k) else as.integer(heads)
  Mbar <- prep$Mbar
  S <- matrix(0, prep$p, prep$p)
  W <- matrix(0, length(hk), prep$H + !is.null(sink))
  for (i in seq_along(hk)) {
    s <- .att_scores(R0, prep, hk[i], score)
    if (!is.null(sink)) s <- c(s, .sink_score(R0, prep, hk[i], score) + sink)
    w <- .normalise(s, tau, normalize)
    if (!is.null(top) && top < length(w)) {
      w[rank(-w, ties.method = "first") > top] <- 0
      w <- w / sum(w)
    }
    Sh <- matrix(prep$V %*% w[seq_len(prep$H)], prep$p, prep$p)
    Sh <- if (sum(diag(Sh)) > 0) Sh / sum(diag(Sh)) else Sh
    if (!is.null(sink)) {
      ws <- w[prep$H + 1L]
      Sh <- (1 - ws) * Sh + ws * Mbar
    }
    S <- S + Sh
    W[i, ] <- w
  }
  S <- S / length(hk)
  if (contrast != 0) S <- S - contrast * Mbar
  if (self != 0) S <- (1 - self) * S + self * R0 / sum(diag(R0))
  list(S = S, W = W)
}

.query_matrix <- function(query) {
  if (inherits(query, "bp_moment")) return(.moment_matrix(query))
  query <- as.matrix(query)
  if (nrow(query) == ncol(query) && isSymmetric(unname(query))) return(query)
  .moment_matrix(energy_moment(query, q = ncol(query)))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Attention weights over transmitted moments
#'
#' Scores each transmitted moment against a query moment and turns the scores
#' into weights. The key of record \eqn{h} is the projector \eqn{P_h} onto its
#' leading \code{k} eigenvectors. With \code{score = "energy"} the score is
#' \eqn{s_h = tr(P_h R_0)}, where \eqn{R_0} is the energy-normalised moment of
#' the query; when \eqn{R_0} is the full moment of a set of beats,
#' \eqn{1 - s_h} is their mean per-beat relative error in the basis of record
#' \eqn{h}. With \code{score = "projection"} the score is
#' \eqn{tr(P_h P_0)/k}, symmetric in the two subspaces, where \eqn{P_0}
#' projects onto the leading \code{k} eigenvectors of \eqn{R_0}. With
#' \code{score = "cosine"} it is the cosine between \eqn{R_0} and \eqn{M_h} in
#' the Frobenius inner product. The scores depend on the eigenvectors only
#' through \eqn{P_h} or \eqn{M_h}, and are therefore unchanged by their signs,
#' order or any rotation within the retained subspace.
#'
#' @param query A \code{bp_moment} object, a \eqn{p \times p} moment matrix or
#'   a beat matrix, from which a full-rank energy moment is computed.
#' @param moments A list of \code{bp_moment} objects.
#' @param k Dimension of the keys.
#' @param tau Temperature. \code{Inf} gives equal weights and \code{0} puts
#'   all weight on the highest score, shared among ties.
#' @param top Optional number of largest weights kept, the others set to zero
#'   and the kept ones renormalised. A compressor then transmits \code{top}
#'   indices and weights instead of one weight per record.
#' @param score One of \code{"energy"}, \code{"projection"} or
#'   \code{"cosine"}.
#' @param normalize \code{"softmax"}, or \code{"sparsemax"}, the Euclidean
#'   projection of \eqn{s/\tau} onto the simplex, which sets low scores to zero.
#' @param sink Optional offset of an extra key, the mean of the transmitted
#'   moments, whose score is its own score plus \code{sink}. Its weight, in the
#'   attribute \code{sink}, goes to the population moment.
#' @return A vector of weights over the records, with the scores in the
#'   attribute \code{scores}, the weight of the extra key in \code{sink} and
#'   the effective number of records \eqn{1 / \sum_h w_h^2} in \code{n_eff}.
#' @examples
#' sim <- simulate_beats(H = 5, n = 60, p = 20, q = 3, seed = 4)
#' m <- lapply(sim$X, energy_moment, q = 5)
#' attention_weights(sim$X[[1]][1:20, ], m[-1], k = 3, tau = 0.01)
#' @export
attention_weights <- function(query, moments, k, tau, top = NULL,
                              score = "energy", normalize = "softmax",
                              sink = NULL) {
  R0 <- .query_matrix(query)
  prep <- .att_prep(moments, k)
  a <- .att_core(R0, prep, k, tau, score, normalize, top, sink)
  w <- a$W[1, ]
  ws <- if (is.null(sink)) 0 else w[length(w)]
  w <- w[seq_len(prep$H)]
  structure(w, scores = .att_scores(R0, prep, k, score), sink = ws,
            n_eff = 1 / sum(w^2))
}

#' Attention pooling of transmitted moments
#'
#' Returns the leading eigenvectors of an attended moment
#' \eqn{S = \sum_h w_h M_h}, where \eqn{M_h} are the rank-q moments already
#' transmitted and \eqn{w_h} are the [attention_weights()] of a query. With
#' \code{tau = Inf} and the defaults it coincides with [pool_energy()]; with
#' \code{tau = 0} it returns the basis of the record of highest score. No
#' communication from the records is added; the query is computed where the
#' beats to be compressed are.
#'
#' Three options change the attended moment, each term scaled to unit trace
#' first. \code{contrast} subtracts \eqn{\lambda} times the mean moment, as in
#' differential attention, so that the directions shared by all records weigh
#' less. \code{self} mixes in the query, \eqn{(1-\alpha) S + \alpha R_0}, as in
#' retrieval interpolation; the basis then depends on the record's own beats
#' and must be sent with them. \code{heads} averages the attended moments of
#' several key dimensions.
#'
#' @inheritParams attention_weights
#' @param q_out Number of components returned. Defaults to \code{k}.
#' @param contrast Weight \eqn{\lambda} of the subtracted mean moment.
#' @param self Weight \eqn{\alpha} of the query moment.
#' @param heads Optional vector of key dimensions, one per head.
#' @return A list with the basis \code{T}, the eigenvalues \code{values}, the
#'   \code{weights} over the records (a matrix with one row per head when
#'   there are several), the weight of the extra key \code{sink}, the
#'   temperature \code{tau} and \code{n_eff}.
#' @examples
#' sim <- simulate_beats(H = 5, n = 60, p = 20, q = 3, seed = 4)
#' m <- lapply(sim$X, energy_moment, q = 5)
#' fit <- pool_attention(sim$X[[1]][1:20, ], m[-1], k = 3, tau = 0.01)
#' round(fit$weights, 3)
#' @export
pool_attention <- function(query, moments, k, tau, q_out = k, top = NULL,
                           score = "energy", normalize = "softmax", sink = NULL,
                           contrast = 0, self = 0, heads = NULL) {
  R0 <- .query_matrix(query)
  prep <- .att_prep(moments, c(k, heads))
  a <- .att_core(R0, prep, k, tau, score, normalize, top, sink, contrast,
                 self, heads)
  e <- eigen((a$S + t(a$S)) / 2, symmetric = TRUE)
  j <- seq_len(as.integer(q_out))
  W <- a$W[, seq_len(prep$H), drop = FALSE]
  list(T = sign_fix(e$vectors[, j, drop = FALSE]), values = e$values[j],
       weights = if (nrow(W) == 1L) as.numeric(W) else W,
       sink = if (is.null(sink)) 0 else a$W[, prep$H + 1L],
       tau = tau, n_eff = 1 / rowSums(W^2), q = length(j), p = prep$p,
       method = "attention")
}

#' Self-attention among records with doubly stochastic weights
#'
#' Replaces each transmitted moment by an average of all of them,
#' \eqn{\tilde M_j = \sum_h A_{jh} M_h}, with \eqn{A = D K D},
#' \eqn{K_{jh} = \exp\{tr(P_j P_h)/(k\tau)\}} and \eqn{D} diagonal, found by
#' Sinkhorn iterations so that every row and column of \eqn{A} sums to one.
#' Because \eqn{K} is symmetric, so is \eqn{A}: each record gives as much
#' weight as it receives, and no record dominates the refined set. The
#' receiver can do this with what it holds.
#'
#' @param moments A list of \code{bp_moment} objects.
#' @param k Dimension of the projectors.
#' @param tau Temperature.
#' @param iters Maximum number of Sinkhorn iterations.
#' @return A list of \code{bp_moment} objects of the same rank, with the
#'   matrix \code{A} in the attribute \code{A}.
#' @examples
#' sim <- simulate_beats(H = 5, n = 60, p = 20, q = 3, seed = 6)
#' m <- lapply(sim$X, energy_moment, q = 5)
#' r <- refine_moments(m, k = 3, tau = 0.05)
#' round(colSums(attr(r, "A")), 6)
#' @export
refine_moments <- function(moments, k, tau, iters = 500L) {
  H <- length(moments)
  keys <- lapply(moments, function(m) m$T[, seq_len(min(k, m$q)), drop = FALSE])
  G <- matrix(0, H, H)
  for (i in seq_len(H)) for (j in i:H)
    G[i, j] <- G[j, i] <- sum(crossprod(keys[[i]], keys[[j]])^2) / k
  K <- exp((G - max(G)) / tau)
  u <- rep(1, H)
  for (it in seq_len(iters)) {
    u1 <- sqrt(u / as.numeric(K %*% u))
    done <- max(abs(u1 - u)) < 1e-12 * max(u1)
    u <- u1
    if (done) break
  }
  A <- K * outer(u, u)
  p <- nrow(moments[[1]]$T)
  V <- vapply(moments, function(m) as.numeric(.moment_matrix(m)), numeric(p * p))
  Vt <- V %*% t(A / rowSums(A))
  out <- lapply(seq_len(H), function(j) {
    M <- matrix(Vt[, j], p, p)
    e <- eigen((M + t(M)) / 2, symmetric = TRUE)
    q <- moments[[j]]$q
    structure(list(T = e$vectors[, seq_len(q), drop = FALSE],
                   values = pmax(e$values[seq_len(q)], 0),
                   tail = moments[[j]]$tail, n = moments[[j]]$n,
                   p = nrow(M), q = q), class = "bp_moment")
  })
  structure(out, A = A)
}

#' Choose attention settings from transmitted moments
#'
#' Leave-one-out over the records, using only what they have transmitted:
#' each record in turn is the query, attends to the others, and the retained
#' error \eqn{1 - tr(T^T M_j T)} of its own rank-q moment is averaged over
#' records. Because \eqn{M_j} is truncated, this error is an upper bound on the
#' per-beat error of the record. By default the query of record \eqn{j} is
#' \eqn{M_j} itself, which is less noisy than a query computed from a few
#' beats and favours smaller temperatures; \code{queries} gives queries
#' computed as they will be for a new record, such as moments of the first
#' beats, at the cost of one more transmitted moment per record.
#'
#' @inheritParams attention_weights
#' @param grid A data frame with one row per candidate and columns among
#'   \code{tau}, \code{score}, \code{normalize}, \code{sink}, \code{contrast},
#'   \code{self}, \code{top}, \code{heads} (a comma-separated string) and
#'   \code{refine} (the temperature of [refine_moments()]). \code{NA} in
#'   \code{sink}, \code{top}, \code{heads} or \code{refine} leaves that option
#'   off. Refinement is recomputed without the left-out record.
#' @param exclude Optional list giving, for each record, the indices of the
#'   records it must not attend to, such as other records of the same patient.
#' @param queries Optional list of query matrices, one per record.
#' @return A list with the chosen row \code{best} and the mean error of every
#'   row, \code{errors}.
#' @examples
#' sim <- simulate_beats(H = 6, n = 60, p = 20, q = 3, seed = 5)
#' m <- lapply(sim$X, energy_moment, q = 5)
#' g <- expand.grid(tau = c(0.01, 0.1, Inf), sink = c(NA, 0))
#' select_attention(m, k = 3, grid = g)$best
#' @export
select_attention <- function(moments, k, grid, exclude = NULL, queries = NULL) {
  grid <- as.data.frame(grid, stringsAsFactors = FALSE)
  get <- function(r, nm, default) {
    if (!nm %in% names(grid)) return(default)
    v <- grid[[nm]][r]
    if (is.na(v)) default else v
  }
  heads_of <- function(r) {
    h <- get(r, "heads", NULL)
    if (is.null(h)) NULL else as.integer(strsplit(as.character(h), ",")[[1]])
  }
  ks <- unique(c(k, unlist(lapply(seq_len(nrow(grid)), heads_of))))
  H <- length(moments)
  Ms <- lapply(moments, .moment_matrix)
  prep <- .att_prep(moments, ks)
  ref_taus <- if ("refine" %in% names(grid))
    unique(grid$refine[!is.na(grid$refine)]) else numeric(0)
  err <- matrix(NA_real_, H, nrow(grid))
  for (j in seq_len(H)) {
    others <- setdiff(seq_len(H), c(j, exclude[[j]]))
    pj <- .att_subset(prep, others)
    pref <- lapply(ref_taus, function(rt)
      .att_prep(refine_moments(moments[others], k, rt), ks))
    names(pref) <- as.character(ref_taus)
    R0 <- if (is.null(queries)) Ms[[j]] else .query_matrix(queries[[j]])
    for (r in seq_len(nrow(grid))) {
      rt <- get(r, "refine", NULL)
      pp <- if (is.null(rt)) pj else pref[[as.character(rt)]]
      a <- .att_core(R0, pp, k, get(r, "tau", Inf), get(r, "score", "energy"),
                     get(r, "normalize", "softmax"), get(r, "top", NULL),
                     get(r, "sink", NULL), get(r, "contrast", 0),
                     get(r, "self", 0), heads_of(r))
      Tj <- .leading(a$S, k)
      err[j, r] <- 1 - sum(Tj * (Ms[[j]] %*% Tj))
    }
  }
  e <- colMeans(err)
  list(best = grid[which.min(e), , drop = FALSE], errors = e)
}

#' Choose the attention temperature from transmitted moments
#'
#' A shortcut for [select_attention()] over temperatures alone.
#'
#' @inheritParams select_attention
#' @param taus Candidate temperatures, \code{Inf} included by default.
#' @return The chosen temperature, with the mean error of each candidate in
#'   the attribute \code{errors}.
#' @examples
#' sim <- simulate_beats(H = 6, n = 60, p = 20, q = 3, seed = 5)
#' m <- lapply(sim$X, energy_moment, q = 5)
#' select_tau(m, k = 3)
#' @export
select_tau <- function(moments, k, taus = c(0, 2^seq(-10, 0), Inf),
                       exclude = NULL, queries = NULL) {
  s <- select_attention(moments, k, data.frame(tau = taus), exclude, queries)
  structure(s$best$tau, errors = stats::setNames(s$errors, taus))
}

#' Attention with a query that grows with the record
#'
#' Splits the beats into consecutive blocks. The basis of each block after the
#' first comes from attention with the energy-normalised moment of all earlier
#' beats as query, so that it uses no beat of the block itself. The
#' transmitted moments play the part of a fixed cache of keys and values, and
#' only the query is updated.
#'
#' @param X A beat matrix, in time order.
#' @param moments A list of \code{bp_moment} objects.
#' @param k Dimension of the bases.
#' @param block Number of beats per block.
#' @param args Arguments of [pool_attention()] other than the query, the
#'   moments and \code{k}, as a list, or a function of the number of earlier
#'   beats returning such a list.
#' @return A list with one basis per block, \code{T}, and the block of each
#'   beat, \code{block}. The first block has no basis.
#' @examples
#' sim <- simulate_beats(H = 5, n = 60, p = 20, q = 3, seed = 7)
#' m <- lapply(sim$X[-1], energy_moment, q = 5)
#' s <- stream_attention(sim$X[[1]], m, k = 3, block = 20,
#'                       args = list(tau = 0.05))
#' length(s$T)
#' @export
stream_attention <- function(X, moments, k, block = 50L, args = list(tau = Inf)) {
  X <- as.matrix(X)
  b <- (seq_len(nrow(X)) - 1L) %/% as.integer(block) + 1L
  nb <- max(b)
  preps <- list()
  Ts <- vector("list", nb)
  for (i in seq_len(nb)[-1L]) {
    past <- X[b < i, , drop = FALSE]
    a <- if (is.function(args)) args(nrow(past)) else args
    key <- paste(c(k, a$heads), collapse = ",")
    if (is.null(preps[[key]])) preps[[key]] <- .att_prep(moments, c(k, a$heads))
    R0 <- .moment_matrix(energy_moment(past, q = ncol(X)))
    s <- .att_core(R0, preps[[key]], k, a$tau %||% Inf, a$score %||% "energy",
                   a$normalize %||% "softmax", a$top, a$sink,
                   a$contrast %||% 0, a$self %||% 0, a$heads)
    Ts[[i]] <- .leading(s$S, k)
  }
  list(T = Ts, block = b)
}
