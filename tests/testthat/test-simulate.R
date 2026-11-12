test_that("simulate_beats returns the documented structure", {
  sim <- simulate_beats(H = 5, n = 80, p = 30, q = 4, epsilon = 0.2,
                        sigma = 0.1, seed = 51)
  expect_length(sim$X, 5)
  expect_true(all(vapply(sim$X, nrow, integer(1)) == 80L))
  expect_true(all(vapply(sim$X, ncol, integer(1)) == 30L))
  expect_length(sim$Q_list, 5)
  expect_equal(dim(sim$Q0), c(30L, 4L))
  expect_equal(length(sim$lambda), 4L)
  expect_equal(sim$params$H, 5L)
  expect_equal(attr(sim$X[[1]], "record"), "sim1")
})

test_that("simulation is reproducible and record sizes can vary", {
  a <- simulate_beats(H = 3, n = c(50, 60, 70), p = 15, q = 3, seed = 52)
  b <- simulate_beats(H = 3, n = c(50, 60, 70), p = 15, q = 3, seed = 52)
  expect_equal(a$X, b$X)
  expect_equal(vapply(a$X, nrow, integer(1)), c(50L, 60L, 70L))
})

test_that("homogeneous records share the common basis", {
  sim <- simulate_beats(H = 4, n = 60, p = 20, q = 3, model = "homogeneous",
                        epsilon = 0, sigma = 0, seed = 53)
  for (h in 1:4) {
    expect_equal(subspace_distance(sim$Q_list[[h]], sim$Q0), 0,
                 tolerance = 1e-6)
  }
})

test_that("heterogeneous records rotate away from the common basis", {
  sim <- simulate_beats(H = 6, n = 60, p = 20, q = 3,
                        model = "heterogeneous", epsilon = 0.5, seed = 54)
  d <- vapply(sim$Q_list, subspace_distance, numeric(1), T2 = sim$Q0)
  expect_true(all(d > 0))
})

test_that("mixture model returns labels and normal indicators", {
  sim <- simulate_beats(H = 2, n = 100, p = 20, q = 3, model = "mixture",
                        mu = 1, pi = 0.4, sigma = 0.1, seed = 55)
  lab <- unlist(sim$labels)
  expect_setequal(unique(lab), c("normal", "abnormal"))
  nrm <- unlist(sim$normal)
  expect_equal(mean(nrm), 0.6, tolerance = 0.15)
  expect_true(all(unlist(sim$labels)[unlist(sim$normal)] == "normal"))
})

test_that("custom spectra and bases are honoured", {
  Q <- random_basis(10, 3, seed = 56)
  sim <- simulate_beats(H = 2, n = 40, p = 10, q = 3, spectrum = c(4, 2, 1),
                        basis = Q, sigma = 0, seed = 57)
  expect_equal(sim$lambda, c(4, 2, 1))
  expect_equal(subspace_distance(sim$Q0, Q), 0, tolerance = 1e-6)
})

test_that("corrupt_basis preserves the span and is reproducible", {
  Q <- random_basis(9, 4, seed = 58)
  for (mode in c("sign", "permutation", "both")) {
    a <- corrupt_basis(Q, mode, seed = 3)
    b <- corrupt_basis(Q, mode, seed = 3)
    expect_equal(a, b)
    expect_equal(subspace_distance(Q, a), 0, tolerance = 1e-6)
  }
})
