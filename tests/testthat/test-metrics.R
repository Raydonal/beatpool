test_that("reconstruction_error is zero for a full basis and monotone in k", {
  set.seed(31)
  X <- matrix(rnorm(60 * 6), 60, 6)
  T <- pca_basis(X, q = 6)$T
  expect_equal(reconstruction_error(X, T), 0, tolerance = 1e-12)
  e <- vapply(1:6, function(k) reconstruction_error(X, T, k), numeric(1))
  expect_true(all(diff(e) <= 1e-12))
  eb <- reconstruction_error(X, T, 3, per_beat = TRUE)
  expect_length(eb, 60)
  expect_equal(mean(eb), reconstruction_error(X, T, 3))
})

test_that("reconstruction_curve returns one row per component", {
  set.seed(32)
  X <- matrix(rnorm(40 * 5), 40, 5)
  T <- random_basis(5, 5, seed = 1)
  cv <- reconstruction_curve(X, T, k = 1:5)
  expect_named(cv, c("k", "error"))
  expect_equal(nrow(cv), 5L)
  expect_true(all(diff(cv$error) <= 1e-12))
})

test_that("compress_beats reconstructs and reports errors", {
  set.seed(33)
  X <- matrix(rnorm(50 * 8), 50, 8)
  T <- pca_basis(X, q = 4)$T
  z <- compress_beats(X, T, k = 4)
  expect_equal(dim(z$Y), c(50L, 4L))
  expect_equal(dim(z$X_hat), dim(X))
  expect_equal(z$mean_error, reconstruction_error(X, T, 4), tolerance = 1e-12)
  expect_equal(z$X_hat, (sweep(X, 2, colMeans(X)) %*% z$T) %*% t(z$T) +
                 rep(colMeans(X), each = nrow(X)), tolerance = 1e-10)
})

test_that("rate_distortion returns a monotone rate--distortion curve", {
  set.seed(37)
  X <- matrix(rnorm(80 * 10), 80, 10)
  T <- pca_basis(X, q = 4)$T
  rd <- rate_distortion(X, T, k = 4, steps = 2^seq(0, -6, by = -1))
  expect_named(rd, c("step", "rate", "distortion", "compression_ratio"))
  expect_equal(nrow(rd), 7L)
  expect_true(all(rd$rate > 0))
  expect_true(all(diff(rd$distortion) < 0))
  expect_equal(rd$compression_ratio, 11 * 10 / rd$rate)
  floor <- 100 * reconstruction_error(X, T, k = 4)
  expect_lt(abs(rd$distortion[nrow(rd)] - floor), 1)
})

test_that("rate_distortion reaches zero distortion on an exact low-rank matrix", {
  set.seed(39)
  Q <- qr.Q(qr(matrix(rnorm(12 * 3), 12, 3)))
  X <- matrix(rnorm(70 * 3), 70, 3) %*% t(Q)
  rd <- rate_distortion(X, Q, k = 3, steps = c(0.5, 0.01))
  expect_lt(rd$distortion[2], 0.05)
  expect_gt(rd$distortion[1], rd$distortion[2])
})

test_that("rate_distortion charges the basis only when it is transmitted", {
  set.seed(38)
  X <- matrix(rnorm(60 * 8), 60, 8)
  T <- random_basis(8, 8, seed = 2)
  a <- rate_distortion(X, T, k = 3, steps = 0.25)
  b <- rate_distortion(X, T, k = 3, steps = 0.25, transmit_basis = FALSE)
  expect_equal(a$rate - b$rate, 3 * 8 * 16 / 60)
  expect_equal(a$distortion, b$distortion)
  c <- rate_distortion(X, T, k = 3, steps = 0.25, n_amortize = 600)
  expect_lt(c$rate, a$rate)
})

test_that("adjusted Rand index matches known values", {
  expect_equal(ari(c(1, 1, 2, 2), c(1, 1, 2, 2)), 1)
  expect_equal(ari(c(1, 1, 2, 2), c(2, 2, 1, 1)), 1)
  expect_equal(ari(c(1, 2, 1, 2), c(1, 2, 2, 1)), -0.5)
  expect_true(is.na(ari(c(1, 1), c(1, 1))))
  expect_true(is.na(ari(c(1, NA, 2), c(1, 2, NA))))
})

test_that("purity and prop_correct agree with hand computations", {
  cl <- c(1, 1, 2, 2)
  tr <- c("N", "N", "N", "V")
  expect_equal(purity(cl, tr), 0.75)
  pc <- prop_correct(rep(1, 5), c("N", "N", "N", "N", "V"))
  expect_equal(pc, 0.8)
  expect_error(prop_correct(1:3, 1:2), "equal length")
})

test_that("cluster_metrics and cluster_beats run", {
  sim <- simulate_beats(H = 1, n = 200, p = 20, q = 3, model = "mixture",
                        mu = 1.5, sigma = 0.1, seed = 34)
  X <- sim$X[[1]]
  T <- pca_basis(X, q = 5)$T
  fit <- cluster_beats(X, T, k = 3, K = 2, truth = sim$labels[[1]], seed = 1)
  expect_length(fit$cluster, 200)
  expect_named(fit$metrics, c("prop_correct", "ari", "purity"))
  expect_gte(fit$metrics$ari, 0.5)
  expect_equal(fit$cluster, cluster_beats(X, T, k = 3, K = 2, seed = 1)$cluster)
})

test_that("bootstrap helpers return coherent intervals", {
  set.seed(35)
  x <- rnorm(60, mean = 2)
  ci <- bootstrap_ci(x, B = 400, seed = 1)
  expect_named(ci, c("lower", "estimate", "upper"))
  expect_lte(ci["lower"], ci["estimate"])
  expect_lte(ci["estimate"], ci["upper"])
  expect_equal(unname(ci["estimate"]), mean(x))
  pb <- paired_bootstrap(x, x, B = 200, seed = 1)
  expect_equal(pb$estimate, 0)
  expect_equal(pb$p_value, 1)
  y <- x + 1
  pb2 <- paired_bootstrap(y, x, B = 400, seed = 1)
  expect_equal(pb2$estimate, 1)
  expect_lt(pb2$p_value, 0.01)
})

test_that("compare_methods returns one row per method", {
  set.seed(36)
  d <- data.frame(a = rnorm(20), b = rnorm(20, 0.2), c = rnorm(20, -0.1))
  out <- compare_methods(d, reference = "a")
  expect_equal(nrow(out), 3L)
  expect_named(out, c("method", "mean", "lower", "upper", "diff",
                      "diff_lower", "diff_upper", "p_value"))
  expect_equal(out$diff[out$method == "a"], 0)
})
