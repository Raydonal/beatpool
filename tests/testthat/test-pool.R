test_that("pooling identical bases recovers the shared subspace", {
  Q <- random_basis(12, 4, seed = 41)
  for (method in c("procrustes", "bar", "qr", "mean")) {
    p <- pool_bases(rep(list(Q), 5), method = method, align = TRUE)
    expect_equal(subspace_distance(p$T, Q), 0, tolerance = 1e-6)
    expect_equal(crossprod(p$T), diag(4), tolerance = 1e-9)
  }
})

test_that("pooling is invariant to sign and permutation when aligned", {
  Q <- random_basis(12, 4, seed = 42)
  bases <- lapply(1:6, function(h)
    corrupt_basis(Q, mode = if (h %% 2) "sign" else "both", seed = h))
  for (method in c("procrustes", "bar", "qr", "mean")) {
    p <- pool_bases(bases, method = method, align = TRUE)
    expect_equal(subspace_distance(p$T, Q), 0, tolerance = 1e-6)
  }
})

test_that("alignment matters when records have heterogeneous bases", {
  sim <- simulate_beats(H = 8, n = 150, p = 20, q = 3, epsilon = 0.4,
                        sigma = 0.05, seed = 43)
  bases <- extract_bases(sim$X, q = 3)
  raw <- pool_bases(bases, method = "bar", align = FALSE)
  aligned <- pool_bases(bases, method = "bar", align = TRUE)
  expect_gt(subspace_distance(raw$T, aligned$T), 1e-3)
  expect_lt(subspace_distance(aligned$T, sim$Q0),
            subspace_distance(raw$T, sim$Q0))
})

test_that("GPA converges and honours the iteration budget", {
  Q <- random_basis(10, 3, seed = 8)
  bases <- lapply(1:8, function(h)
    orthonormalize(Q + matrix(rnorm(30, sd = 0.05), 10, 3), 3))
  p <- pool_bases(bases, iterations = 10)
  expect_true(isTRUE(p$converged))
  expect_lte(p$iterations_used, 10L)
  g <- gpa_bases(bases)
  expect_s3_class(g, "bp_pool")
  expect_lt(subspace_distance(g$T, p$T), 1e-6)
})

test_that("weights are normalised and replicate weighting", {
  Qs <- lapply(1:4, function(h) random_basis(8, 3, seed = h))
  p1 <- pool_bases(Qs, weights = c(1, 1, 1, 1))
  p2 <- pool_bases(Qs, weights = c(2, 2, 2, 2))
  expect_equal(p1$T, p2$T, tolerance = 1e-12)
  expect_equal(sum(p1$weights), 1)
  expect_error(pool_bases(Qs, weights = 1:2), "wrong length")
})

test_that("q truncates the pooled basis", {
  Qs <- lapply(1:5, function(h) random_basis(9, 5, seed = h))
  p <- pool_bases(Qs, q = 2)
  expect_equal(ncol(p$T), 2L)
  expect_equal(p$q, 2L)
})

test_that("extract_bases and pool_from_beats work end to end", {
  sim <- simulate_beats(H = 4, n = 120, p = 20, q = 4, epsilon = 0.1,
                        sigma = 0.1, seed = 9)
  bases <- extract_bases(sim$X, q = 4)
  expect_length(bases, 4)
  expect_true(all(vapply(bases, inherits, logical(1), "bp_basis")))
  p <- pool_from_beats(sim$X, q = 4)
  expect_s3_class(p, "bp_pool")
  expect_equal(ncol(p$T), 4L)
  expect_equal(length(p$values), 4L)
  expect_true(all(diff(p$values) <= 1e-8))
})

test_that("leave-one-out pooling and jackknife stability", {
  sim <- simulate_beats(H = 5, n = 90, p = 15, q = 3, epsilon = 0.2,
                        sigma = 0.1, seed = 10)
  bases <- extract_bases(sim$X, q = 3)
  loo <- leave_one_out_pool(bases, data = sim$X)
  expect_length(loo$leave_one_out, 5)
  st <- jackknife_stability(loo)
  expect_equal(nrow(st), 5L)
  expect_named(st, c("record", "relative", "subspace"))
  expect_true(all(st$subspace >= 0))
})

test_that("robust pooling downweights an outlying record", {
  Q <- random_basis(10, 3, seed = 12)
  bases <- lapply(1:6, function(h)
    orthonormalize(Q + matrix(rnorm(30, sd = 0.05), 10, 3), 3))
  bases[[6]] <- random_basis(10, 3, seed = 999)
  p <- pool_bases(bases, robust = TRUE)
  expect_equal(sum(p$weights), 1, tolerance = 1e-8)
  expect_lt(p$weights[6], min(p$weights[-6]))
})

test_that("hierarchical pooling returns one basis per dimension", {
  sim <- simulate_beats(H = 5, n = 100, p = 18, q = 4, epsilon = 0.15,
                        sigma = 0.1, seed = 13)
  bases <- extract_bases(sim$X, q = 4)
  hp <- hierarchical_pool(bases, q = 3)
  expect_s3_class(hp, "bp_hpool")
  expect_length(hp, 3)
  expect_equal(unname(vapply(hp, function(z) ncol(z$T), integer(1))), 1:3)
  for (k in 1:3) {
    expect_lt(subspace_distance(hp[[k]]$T, sim$Q0[, 1:k, drop = FALSE]), 0.5)
  }
})

test_that("pool_spectral is exact under sign and permutation corruption", {
  sim <- simulate_beats(H = 6, n = 150, p = 25, q = 5, epsilon = 0.1,
                        sigma = 0.05, seed = 14)
  bases <- extract_bases(sim$X, q = 5)
  corrupted <- lapply(seq_along(bases), function(h) {
    b <- bases[[h]]
    p <- sample.int(5)
    s <- sample(c(-1, 1), 5, replace = TRUE)
    b$T <- sweep(b$T[, p, drop = FALSE], 2L, s, "*")
    b$values <- b$values[p]
    b
  })
  s1 <- pool_spectral(bases, q = 5, q_out = 5, return_cov = TRUE)
  s2 <- pool_spectral(corrupted, q = 5, q_out = 5, return_cov = TRUE)
  expect_equal(s1$S, s2$S, tolerance = 1e-10)
  expect_equal(subspace_distance(s1$T, s2$T), 0, tolerance = 1e-6)
})

test_that("pool_spectral reproduces covariance pooling for full bases", {
  sim <- simulate_beats(H = 5, n = 120, p = 12, q = 4, epsilon = 0.2,
                        sigma = 0.05, seed = 15)
  X <- sim$X
  p <- ncol(X[[1]])
  bases <- extract_bases(X, q = p)
  sp <- pool_spectral(bases, q = p, q_out = p, return_cov = TRUE)
  cv <- pool_covariance(X, q = p, return_cov = TRUE)
  expect_equal(sp$S, cv$S, tolerance = 1e-8)
  expect_equal(sp$values, cv$values, tolerance = 1e-8)
  expect_lt(subspace_distance(sp$T, cv$T), 1e-8)
})

test_that("residual imputation preserves the trace", {
  sim <- simulate_beats(H = 4, n = 100, p = 10, q = 3, epsilon = 0.1,
                        sigma = 0.05, seed = 16)
  bases <- extract_bases(sim$X, q = 3)
  S <- pool_spectral(bases, q = 3, q_out = 3, residual = TRUE,
                     return_cov = TRUE)$S
  target <- mean(vapply(bases, function(b) b$total, numeric(1)))
  expect_equal(sum(diag(S)), target, tolerance = 1e-8)
})

test_that("print methods return their input invisibly", {
  sim <- simulate_beats(H = 3, n = 80, p = 12, q = 3, seed = 17)
  b <- pca_basis(sim$X[[1]], q = 3)
  expect_invisible(print(b))
  p <- pool_bases(list(b, b))
  expect_invisible(print(p))
  hp <- hierarchical_pool(list(b, b), q = 2)
  expect_invisible(print(hp))
  expect_equal(as.matrix(p), p$T)
})
