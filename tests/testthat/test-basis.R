test_that("sign_fix makes the largest absolute entry positive", {
  T <- random_basis(10, 2, seed = 1)
  F <- sign_fix(T)
  for (j in 1:2) {
    v <- F[, j]
    expect_gte(v[which.max(abs(v))], 0)
    expect_equal(abs(sum(T[, j] * v)), 1, tolerance = 1e-12)
  }
})

test_that("pca_basis returns orthonormal columns and accounted variance", {
  set.seed(2)
  X <- matrix(rnorm(100 * 6), 100, 6)
  b <- pca_basis(X, q = 3)
  expect_s3_class(b, "bp_basis")
  expect_equal(crossprod(b$T), diag(3), tolerance = 1e-10)
  expect_equal(ncol(b$T), 3L)
  expect_equal(length(b$values), 3L)
  expect_lte(sum(b$values), b$total + 1e-10)
  expect_true(all(diff(b$values) <= 1e-12))
})

test_that("pca_basis with q = p reconstructs exactly", {
  set.seed(3)
  X <- matrix(rnorm(50 * 4), 50, 4)
  b <- pca_basis(X, q = 4, scores = TRUE)
  expect_equal(reconstruction_error(X, b$T), 0, tolerance = 1e-10)
  expect_equal(sum(b$values), b$total, tolerance = 1e-10)
  expect_equal(dim(b$scores), c(50L, 4L))
})

test_that("dhmt_basis is orthogonal and square", {
  T <- dhmt_basis(15)
  expect_equal(dim(T), c(15L, 15L))
  expect_equal(crossprod(T), diag(15), tolerance = 1e-10)
  expect_equal(T, dhmt_basis(15))
})

test_that("dct_basis is orthogonal with a constant first column", {
  T <- dct_basis(9)
  expect_equal(crossprod(T), diag(9), tolerance = 1e-10)
  expect_equal(T[, 1], rep(1 / 3, 9), tolerance = 1e-12)
})

test_that("random_basis has orthonormal columns and honours the seed", {
  A <- random_basis(8, 3, seed = 7)
  B <- random_basis(8, 3, seed = 7)
  expect_equal(A, B)
  expect_equal(crossprod(A), diag(3), tolerance = 1e-12)
})

test_that("align_basis removes sign and permutation indeterminacy", {
  Q <- random_basis(10, 4, seed = 11)
  for (mode in c("sign", "permutation", "both")) {
    Qc <- corrupt_basis(Q, mode = mode, seed = 5)
    al <- align_basis(Qc, Q)
    expect_equal(subspace_distance(al$Q, Q), 0, tolerance = 1e-6)
    expect_equal(sum(sin(principal_angles(al$Q, Q))^2), 0, tolerance = 1e-10)
  }
})

test_that("principal angles and subspace distances have the right scale", {
  Q <- random_basis(12, 3, seed = 21)
  expect_equal(principal_angles(Q, Q), rep(0, 3), tolerance = 1e-6)
  expect_equal(subspace_distance(Q, Q), 0, tolerance = 1e-6)
  expect_equal(subspace_distance(Q, Q, "max"), 0, tolerance = 1e-6)
  expect_equal(subspace_distance(Q, Q, "frobenius"), 0, tolerance = 1e-6)
  Qc <- corrupt_basis(Q, "both", seed = 2)
  expect_equal(subspace_distance(Q, Qc), 0, tolerance = 1e-6)
  O <- random_basis(12, 3, seed = 33)
  expect_gt(subspace_distance(Q, O), 0)
  expect_equal(subspace_distance(Q, O), subspace_distance(O, Q))
})

test_that("orthonormalize and complete_basis span the same space", {
  set.seed(4)
  A <- matrix(rnorm(30), 10, 3)
  Q <- orthonormalize(A)
  expect_equal(ncol(Q), 3L)
  expect_equal(crossprod(Q), diag(3), tolerance = 1e-12)
  expect_equal(subspace_distance(Q, A), 0, tolerance = 1e-6)
  Z <- complete_basis(Q)
  expect_equal(dim(Z), c(10L, 10L))
  expect_equal(crossprod(Z), diag(10), tolerance = 1e-10)
  expect_equal(subspace_distance(Z[, 1:3], Q), 0, tolerance = 1e-6)
})

test_that("basis_variance and order_basis agree with direct computation", {
  set.seed(5)
  X <- matrix(rnorm(80 * 5), 80, 5)
  b <- pca_basis(X, q = 5)
  v <- basis_variance(b$T, list(X))
  expect_equal(v, b$values, tolerance = 1e-8)
  o <- order_basis(b$T[, c(4, 2, 5, 1, 3)], list(X))
  expect_equal(o$order, c(4, 2, 5, 1, 3))
  expect_equal(subspace_distance(o$T, b$T), 0, tolerance = 1e-6)
})

test_that("basis_frame returns long format", {
  T <- random_basis(6, 2, seed = 1)
  f <- basis_frame(T)
  expect_equal(nrow(f), 12L)
  expect_named(f, c("time", "component", "value"))
  expect_equal(f$value, as.vector(T))
})
