test_that("the energy moment reproduces the per-beat error", {
  sim <- simulate_beats(H = 1, n = 80, p = 16, q = 3, seed = 51)
  X <- sim$X[[1]]
  m <- energy_moment(X, q = 16)
  Q <- random_basis(16, 4, seed = 52)
  R <- m$T %*% (m$values * t(m$T))
  expect_equal(reconstruction_error(X, Q), 1 - sum(diag(t(Q) %*% R %*% Q)),
               tolerance = 1e-10)
  expect_equal(sum(m$values) + m$tail, 1, tolerance = 1e-10)
})

test_that("the pooled energy moment ignores sign and permutation", {
  sim <- simulate_beats(H = 5, n = 60, p = 16, q = 3, seed = 53)
  m <- lapply(sim$X, energy_moment, q = 5)
  a <- pool_energy(m, q_out = 3)
  m2 <- lapply(m, function(x) {
    o <- rev(seq_len(x$q))
    x$T <- -x$T[, o]
    x$values <- x$values[o]
    x
  })
  expect_equal(subspace_distance(a$T, pool_energy(m2, q_out = 3)$T), 0,
               tolerance = 1e-8)
})

test_that("one shared basis is the pooled energy basis", {
  sim <- simulate_beats(H = 6, n = 60, p = 16, q = 3, seed = 54)
  m <- lapply(sim$X, energy_moment, q = 5)
  fit <- pool_shared(m, k = 3, C = 1, starts = 1)
  expect_equal(subspace_distance(fit$T[[1]], pool_energy(m, q_out = 3)$T), 0,
               tolerance = 1e-8)
})

test_that("the objective does not fall with more bases", {
  sim <- simulate_beats(H = 9, n = 60, p = 16, q = 3, epsilon = 0.5, seed = 55)
  m <- lapply(sim$X, energy_moment, q = 5)
  v <- vapply(1:3, function(C) pool_shared(m, k = 2, C = C, seed = 1)$value,
              numeric(1))
  expect_true(all(diff(v) >= -1e-10))
})

test_that("separated groups are recovered and selected", {
  Qs <- list(diag(12)[, 1:2], diag(12)[, 7:8])
  gen <- function(Q, seed) {
    set.seed(seed)
    matrix(rnorm(400), 200, 2) %*% t(Q) + matrix(rnorm(2400, sd = 0.02), 200, 12)
  }
  X <- lapply(1:8, function(h) gen(Qs[[1 + h %% 2]], h))
  fit <- pool_shared(lapply(X, energy_moment, q = 3), k = 2, C = 2, seed = 1)
  truth <- 1 + (1:8) %% 2
  expect_equal(length(unique(paste(fit$group, truth))), 2L)
  held <- gen(Qs[[1]], 99)
  i <- select_shared(fit, held)
  expect_lt(subspace_distance(fit$T[[i]], Qs[[1]]), 0.1)
})

test_that("a seed leaves the random state alone", {
  sim <- simulate_beats(H = 4, n = 40, p = 12, q = 2, seed = 56)
  m <- lapply(sim$X, energy_moment, q = 3)
  set.seed(1); a <- runif(1)
  set.seed(1); invisible(pool_shared(m, k = 2, C = 2, starts = 1, seed = 9)); b <- runif(1)
  expect_identical(a, b)
  expect_error(pool_shared(m, k = 2, C = 9), "cannot exceed")
})
