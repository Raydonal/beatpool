test_that("segment_beats extracts windows of the right length", {
  set.seed(41)
  x <- as.numeric(arima.sim(list(ar = 0.7), n = 600))
  peaks <- c(100, 200, 300, 400, 500)
  X <- segment_beats(x, peaks, wsize = 10, normalize = "none")
  expect_equal(dim(X), c(5L, 21L))
  expect_equal(X[1, ], x[90:110])
  expect_equal(attr(X, "peaks"), peaks)
  expect_equal(attr(X, "wsize"), 10L)
  expect_equal(attr(X, "normalize"), "none")
})

test_that("segment_beats handles labels, sorting and duplicates", {
  x <- rep(sin(seq_len(400) / 10), 2)
  peaks <- c(300, 100, 200, 100, 250)
  X <- segment_beats(x, peaks, labels = c("a", "b", "c", "d", "e"),
                     normal = c(FALSE, TRUE, TRUE, TRUE, FALSE),
                     wsize = 5)
  expect_equal(nrow(X), 4L)
  expect_equal(attr(X, "peaks"), c(100, 200, 250, 300))
  expect_equal(attr(X, "labels"), c("b", "c", "e", "a"))
  expect_equal(attr(X, "normal"), c(TRUE, TRUE, FALSE, FALSE))
})

test_that("segment_beats normalisation produces unit-scale windows", {
  set.seed(42)
  x <- cumsum(rnorm(500))
  peaks <- seq(100, 400, by = 50)
  X <- segment_beats(x, peaks, wsize = 12, normalize = "qrs")
  expect_equal(unname(rowMeans(X)), rep(0, nrow(X)), tolerance = 1e-10)
  expect_equal(apply(X, 1, sd), rep(1, nrow(X)), tolerance = 1e-10)
})

test_that("drop_first_last and edges are respected", {
  x <- rnorm(300)
  peaks <- c(2, 50, 100, 150, 298)
  X <- segment_beats(x, peaks, wsize = 10, drop_edges = TRUE,
                     drop_first_last = TRUE)
  expect_equal(attr(X, "peaks"), c(50, 100, 150))
  expect_equal(nrow(segment_beats(x, peaks, wsize = 10)), 3L)
  expect_error(segment_beats(x, peaks, wsize = 0, normalize = "none"), NA)
})

test_that("empty input gives an empty matrix with attributes", {
  X <- segment_beats(rep(0, 100), integer(0))
  expect_equal(dim(X), c(0L, 61L))
  expect_equal(attr(X, "peaks"), integer(0))
})

test_that("prepare_record keeps only annotated beats", {
  ann <- data.frame(
    index = c(50, 100, 150, 200),
    symbol = c("N", "N", "V", "N"),
    beat = c(TRUE, TRUE, FALSE, TRUE),
    normal = c(TRUE, TRUE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  x <- sin(seq_len(400) / 7)
  X <- prepare_record(x, ann, record = "test", fs = 360, wsize = 10,
                      drop_first_last = FALSE)
  expect_equal(nrow(X), 3L)
  expect_equal(attr(X, "labels"), c("N", "N", "N"))
  expect_equal(attr(X, "record"), "test")
  expect_equal(attr(X, "lead"), "MLII")
})

test_that("beat_summary reports one row per beat", {
  X <- segment_beats(rnorm(300), c(100, 150, 200), wsize = 10)
  s <- beat_summary(X)
  expect_equal(nrow(s), 3L)
  expect_named(s, c("record", "beat", "label", "normal", "mean", "sd",
                    "min", "max", "energy"))
  expect_equal(s$beat, 1:3)
})

test_that("read_mitdb_csv selects the requested lead", {
  f <- tempfile(fileext = ".csv")
  writeLines(c("MIT-BIH export",
               "time,MLII,V5",
               "0.000,1.0,2.0",
               "0.003,1.5,2.5"), f)
  d1 <- read_mitdb_csv(f, lead = "MLII")
  d2 <- read_mitdb_csv(f, lead = "V5")
  expect_equal(d1$x, c(1.0, 1.5))
  expect_equal(d2$x, c(2.0, 2.5))
  expect_equal(d1$time, c(0, 0.003))
  unlink(f)
})
