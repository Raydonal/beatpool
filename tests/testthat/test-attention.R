test_that("the score is one minus the per-beat error in the record basis", {
  sim <- simulate_beats(H = 4, n = 80, p = 16, q = 3, seed = 61)
  m <- lapply(sim$X, energy_moment, q = 5)
  X0 <- sim$X[[1]]
  w <- attention_weights(X0, m[-1], k = 3, tau = 0.1)
  e <- vapply(m[-1], function(z) reconstruction_error(X0, z$T[, 1:3], k = 3),
              numeric(1))
  expect_equal(attr(w, "scores"), 1 - e, tolerance = 1e-10)
  expect_equal(sum(w), 1)
})

test_that("infinite temperature gives the pooled energy moment", {
  sim <- simulate_beats(H = 5, n = 60, p = 16, q = 3, seed = 62)
  m <- lapply(sim$X, energy_moment, q = 5)
  a <- pool_attention(sim$X[[1]][1:20, ], m, k = 3, tau = Inf)
  expect_equal(a$weights, rep(0.2, 5))
  expect_equal(subspace_distance(a$T, pool_energy(m, q_out = 3)$T), 0,
               tolerance = 1e-6)
})

test_that("zero temperature returns the best record", {
  sim <- simulate_beats(H = 5, n = 60, p = 16, q = 3, seed = 63)
  m <- lapply(sim$X, energy_moment, q = 5)
  a <- pool_attention(sim$X[[2]], m, k = 3, tau = 0)
  expect_equal(which.max(a$weights), 2L)
  expect_equal(a$n_eff, 1)
  expect_equal(subspace_distance(a$T, m[[2]]$T[, 1:3]), 0, tolerance = 1e-6)
})

test_that("attention ignores sign, order and rotation of the eigenvectors", {
  sim <- simulate_beats(H = 5, n = 60, p = 16, q = 3, seed = 64)
  m <- lapply(sim$X, energy_moment, q = 3)
  q0 <- sim$X[[1]][1:30, ]
  a <- pool_attention(q0, m, k = 3, tau = 0.05)
  rot <- random_basis(3, 3, seed = 65)
  m2 <- lapply(m, function(x) {
    o <- rev(seq_len(x$q))
    x$T <- -x$T[, o]
    x$values <- x$values[o]
    x
  })
  b <- pool_attention(q0, m2, k = 3, tau = 0.05)
  expect_equal(a$weights, b$weights, tolerance = 1e-10)
  expect_equal(subspace_distance(a$T, b$T), 0, tolerance = 1e-6)
  w1 <- attention_weights(q0, m, k = 3, tau = 0.05)
  m3 <- lapply(m, function(x) { x$T <- x$T %*% rot; x })
  w3 <- attention_weights(q0, m3, k = 3, tau = 0.05)
  expect_equal(as.numeric(w1), as.numeric(w3), tolerance = 1e-10)
})

test_that("the temperature is chosen from the candidates", {
  sim <- simulate_beats(H = 6, n = 60, p = 16, q = 3, seed = 66)
  m <- lapply(sim$X, energy_moment, q = 5)
  taus <- c(0, 0.01, Inf)
  t0 <- select_tau(m, k = 3, taus = taus)
  expect_true(t0 %in% taus)
  expect_length(attr(t0, "errors"), 3L)
})

test_that("sparse attention keeps the largest weights", {
  sim <- simulate_beats(H = 6, n = 60, p = 16, q = 3, seed = 67)
  m <- lapply(sim$X, energy_moment, q = 5)
  w <- attention_weights(sim$X[[1]], m, k = 3, tau = 0.05)
  w2 <- attention_weights(sim$X[[1]], m, k = 3, tau = 0.05, top = 2)
  expect_equal(sum(w2 > 0), 2L)
  expect_equal(sum(w2), 1)
  expect_equal(which(w2 > 0), sort(order(-w)[1:2]))
})

test_that("the temperature search matches direct attention pooling", {
  sim <- simulate_beats(H = 5, n = 60, p = 12, q = 3, seed = 68)
  m <- lapply(sim$X, energy_moment, q = 4)
  taus <- c(0, 0.02, 0.2, Inf)
  e <- attr(select_tau(m, k = 2, taus = taus), "errors")
  direct <- vapply(taus, function(tau) mean(vapply(seq_along(m), function(j) {
    Mj <- m[[j]]$T %*% (m[[j]]$values * t(m[[j]]$T))
    Tj <- pool_attention(Mj, m[-j], k = 2, tau = tau)$T
    1 - sum(Tj * (Mj %*% Tj))
  }, numeric(1))), numeric(1))
  expect_equal(unname(e), direct, tolerance = 1e-8)
})

test_that("sparsemax weights lie on the simplex and can be exactly zero", {
  sim <- simulate_beats(H = 8, n = 60, p = 16, q = 3, model = "heterogeneous",
                        seed = 69)
  m <- lapply(sim$X, energy_moment, q = 5)
  w <- attention_weights(sim$X[[1]], m, k = 3, tau = 0.01,
                         normalize = "sparsemax")
  expect_equal(sum(w), 1)
  expect_true(all(w >= 0))
  expect_true(any(w == 0))
})

test_that("a dominant sink returns the pooled energy basis", {
  sim <- simulate_beats(H = 5, n = 60, p = 16, q = 3, seed = 70)
  m <- lapply(sim$X, energy_moment, q = 5)
  a <- pool_attention(sim$X[[1]], m, k = 3, tau = 0.01, sink = 100)
  expect_equal(a$sink, 1)
  expect_equal(subspace_distance(a$T, pool_energy(m, q_out = 3)$T), 0,
               tolerance = 1e-6)
})

test_that("full self weight returns the basis of the query", {
  sim <- simulate_beats(H = 5, n = 60, p = 16, q = 3, seed = 71)
  m <- lapply(sim$X[-1], energy_moment, q = 5)
  a <- pool_attention(sim$X[[1]], m, k = 3, tau = 0.05, self = 1)
  own <- energy_moment(sim$X[[1]], q = 3)$T
  expect_equal(subspace_distance(a$T, own), 0, tolerance = 1e-6)
})

test_that("one head equals no heads and the symmetric score is symmetric", {
  sim <- simulate_beats(H = 5, n = 60, p = 16, q = 3, seed = 72)
  m <- lapply(sim$X, energy_moment, q = 5)
  q0 <- sim$X[[1]][1:30, ]
  a <- pool_attention(q0, m, k = 3, tau = 0.05)
  b <- pool_attention(q0, m, k = 3, tau = 0.05, heads = 3)
  expect_equal(subspace_distance(a$T, b$T), 0, tolerance = 1e-6)
  s12 <- attr(attention_weights(m[[1]], m[2], k = 3, tau = 1,
                                score = "projection"), "scores")
  s21 <- attr(attention_weights(m[[2]], m[1], k = 3, tau = 1,
                                score = "projection"), "scores")
  expect_equal(s12, s21, tolerance = 1e-10)
})

test_that("Sinkhorn self-attention is symmetric and doubly stochastic", {
  sim <- simulate_beats(H = 6, n = 60, p = 16, q = 3, seed = 73)
  m <- lapply(sim$X, energy_moment, q = 5)
  r <- refine_moments(m, k = 3, tau = 0.05)
  A <- attr(r, "A")
  expect_equal(A, t(A), tolerance = 1e-10)
  expect_equal(rowSums(A), rep(1, 6), tolerance = 1e-8)
  expect_length(r, 6L)
  expect_s3_class(r[[1]], "bp_moment")
})

test_that("streaming bases use only earlier beats", {
  sim <- simulate_beats(H = 5, n = 60, p = 16, q = 3, seed = 74)
  m <- lapply(sim$X[-1], energy_moment, q = 5)
  X <- sim$X[[1]]
  s <- stream_attention(X, m, k = 3, block = 20, args = list(tau = 0.05))
  expect_null(s$T[[1]])
  T2 <- pool_attention(X[1:20, ], m, k = 3, tau = 0.05)$T
  expect_equal(subspace_distance(s$T[[2]], T2), 0, tolerance = 1e-6)
  X2 <- X
  X2[41:60, ] <- X2[41:60, ] * 3 + 1
  s2 <- stream_attention(X2, m, k = 3, block = 20, args = list(tau = 0.05))
  expect_equal(subspace_distance(s$T[[3]], s2$T[[3]]), 0, tolerance = 1e-8)
})

test_that("the settings search returns a row of the grid", {
  sim <- simulate_beats(H = 6, n = 60, p = 16, q = 3, seed = 75)
  m <- lapply(sim$X, energy_moment, q = 5)
  g <- data.frame(tau = c(0.05, 0.05, Inf), sink = c(NA, 0, NA),
                  refine = c(NA, NA, 0.1))
  s <- select_attention(m, k = 3, grid = g)
  expect_equal(nrow(s$best), 1L)
  expect_length(s$errors, 3L)
})

test_that("Sinkhorn self-attention leaves the mean moment unchanged", {
  sim <- simulate_beats(H = 6, n = 60, p = 16, q = 3, seed = 76)
  m <- lapply(sim$X, energy_moment, q = 16)
  r <- refine_moments(m, k = 3, tau = 0.05)
  mean_of <- function(x) Reduce(`+`, lapply(x, function(z)
    z$T %*% (z$values * t(z$T)))) / length(x)
  expect_equal(mean_of(r), mean_of(m), tolerance = 1e-6)
})
