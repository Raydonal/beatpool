test_that("pool_all returns every operator with the right dimensions", {
  sim <- simulate_beats(H = 6, n = 120, p = 20, q = 3, seed = 11)
  all <- pool_all(sim$X, q = 3)
  expect_s3_class(all, "bp_poolall")
  expect_true(all(c("bar", "bara", "qr", "qra", "proc", "gpa", "hier",
                    "spec", "cov", "covn") %in% names(all)))
  expect_true(all(vapply(all, nrow, integer(1)) == 20L))
  expect_true(all(vapply(all, ncol, integer(1)) == 3L))
  num <- attr(all, "numbers")
  expect_equal(unname(num[["bar"]]), 60)
  expect_equal(unname(num[["spec"]]), 3 * 21)
  expect_equal(unname(num[["cov"]]), 20 * 21 / 2)
})

test_that("pool_all without data omits the covariance operators", {
  sim <- simulate_beats(H = 5, n = 100, p = 12, q = 2, seed = 12)
  all <- pool_all(bases = extract_bases(sim$X, q = 2), q = 2)
  expect_false("cov" %in% names(all))
  expect_true("spec" %in% names(all))
})

test_that("summary.bp_pool reports a single row", {
  sim <- simulate_beats(H = 5, n = 100, p = 12, q = 2, seed = 13)
  p <- pool_from_beats(sim$X, q = 2)
  s <- summary(p)
  expect_s3_class(s, "data.frame")
  expect_equal(nrow(s), 1L)
  expect_equal(s$components, 2L)
  expect_true(is.na(s$top_variance) || s$top_variance > 0)
})

test_that("plot_cost_error needs a named error vector", {
  skip_if_not_installed("ggplot2")
  sim <- simulate_beats(H = 5, n = 100, p = 12, q = 2, seed = 14)
  all <- pool_all(sim$X, q = 2)
  err <- vapply(all, function(T) mean(vapply(sim$X, reconstruction_error,
    numeric(1), T = T, k = 2)), numeric(1))
  p <- plot_cost_error(all, err)
  expect_s3_class(p, "ggplot")
  expect_error(plot_cost_error(all, unname(err)), "named")
})

test_that("input validation rejects malformed matrices", {
  expect_error(pca_basis("not a matrix"), "numeric matrix")
  X <- matrix(1:20, 4, 5)
  X[1, 1] <- NA
  expect_error(pca_basis(X), "non-finite")
})
