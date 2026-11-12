word_bytes <- function(words) {
  low <- words %% 256L
  high <- words %/% 256L
  as.raw(as.integer(rbind(low, high)))
}

write_header <- function(dir, record, format, gain = "200", n_sig = 2L,
                         fs = 360, n_samples = 100L,
                         desc = c("MLII", "V5")) {
  first <- paste(record, n_sig, fs, n_samples)
  sigs <- vapply(seq_len(n_sig), function(ch) {
    paste0(record, ".dat", " ", format, " ", gain, " 11 0 0 0 0 ", desc[ch])
  }, character(1))
  writeLines(c(first, sigs), file.path(dir, paste0(record, ".hea")))
}

test_that("ann_codes maps symbols and beat flags", {
  expect_equal(ann_codes$symbol[ann_codes$code == 1], "N")
  expect_equal(ann_codes$symbol[ann_codes$code == 5], "V")
  expect_true(ann_codes$beat[ann_codes$code == 1])
  expect_true(ann_codes$normal[ann_codes$code == 1])
  expect_false(ann_codes$normal[ann_codes$code == 5])
  expect_false(ann_codes$beat[ann_codes$code == 14])
  expect_equal(nrow(ann_codes), 30L)
})

test_that("read_wfdb_header parses the standard layout", {
  d <- tempfile()
  dir.create(d)
  write_header(d, "rec", 16L, gain = "200(0)/mV", n_samples = 50L)
  h <- read_wfdb_header(file.path(d, "rec.hea"))
  expect_equal(h$record, "rec")
  expect_equal(h$n_sig, 2L)
  expect_equal(h$fs, 360)
  expect_equal(h$n_samples, 50L)
  expect_equal(h$signals$gain, c(200, 200))
  expect_equal(h$signals$baseline, c(0, 0))
  expect_equal(h$signals$description, c("MLII", "V5"))
  unlink(d, recursive = TRUE)
})

test_that("read_wfdb decodes format 16 with gain and baseline", {
  d <- tempfile()
  dir.create(d)
  write_header(d, "rec", 16L, n_sig = 2L, n_samples = 10L)
  v1 <- 0:9
  v2 <- 1000:1009
  con <- file(file.path(d, "rec.dat"), "wb")
  writeBin(word_bytes(as.integer(rbind(v1, v2))), con)
  close(con)
  x <- read_wfdb("rec", dir = d)
  expect_equal(nrow(x), 10L)
  expect_equal(names(x), c("time", "MLII", "V5"))
  expect_equal(x$MLII, v1 / 200)
  expect_equal(x$V5, v2 / 200)
  expect_equal(x$time[2], 1 / 360)
  one <- read_wfdb("rec", dir = d, channels = 2L, from = 3, to = 5)
  expect_equal(names(one), c("time", "V5"))
  expect_equal(one$V5, v2[3:5] / 200)
  unlink(d, recursive = TRUE)
})

test_that("read_wfdb treats a zero gain as unity", {
  d <- tempfile()
  dir.create(d)
  write_header(d, "rec", 16L, gain = "0", n_sig = 1L, n_samples = 4L,
               desc = "ECG1")
  con <- file(file.path(d, "rec.dat"), "wb")
  writeBin(word_bytes(c(10L, 20L, 30L, 40L)), con)
  close(con)
  x <- read_wfdb("rec", dir = d)
  expect_equal(x$ECG1, c(10, 20, 30, 40))
  unlink(d, recursive = TRUE)
})

test_that("read_wfdb decodes format 212", {
  d <- tempfile()
  dir.create(d)
  write_header(d, "rec", 212L, n_sig = 1L, fs = 128, n_samples = 2L,
               desc = "ECG1")
  con <- file(file.path(d, "rec.dat"), "wb")
  writeBin(as.raw(c(100, 176, 255)), con)
  close(con)
  x <- read_wfdb("rec", dir = d)
  expect_equal(x$ECG1, c(100, -5) / 200)
  unlink(d, recursive = TRUE)
})

test_that("read_wfdb rejects unsupported formats", {
  d <- tempfile()
  dir.create(d)
  write_header(d, "rec", 24L, n_sig = 1L, n_samples = 4L, desc = "ECG1")
  expect_error(read_wfdb("rec", dir = d), "unsupported")
  unlink(d, recursive = TRUE)
})

test_that("read_annotations handles times, aux strings and skips", {
  d <- tempfile()
  dir.create(d)
  words <- c(1L * 1024L + 10L,
             5L * 1024L + 5L,
             63L * 1024L + 4L, 1L, 1L,
             59L * 1024L, 1000L, 0L,
             1L * 1024L, 0L,
             12L * 1024L + 3L)
  con <- file(file.path(d, "rec.atr"), "wb")
  writeBin(word_bytes(words), con)
  close(con)
  a <- read_annotations(file.path(d, "rec.atr"))
  expect_equal(a$sample, c(10L, 15L, 1000L, 1003L))
  expect_equal(a$code, c(1L, 5L, 1L, 12L))
  expect_equal(a$symbol, c("N", "V", "N", "/"))
  expect_equal(a$normal, c(TRUE, FALSE, TRUE, FALSE))
  expect_equal(a$beat, c(TRUE, TRUE, TRUE, TRUE))
  unlink(d, recursive = TRUE)
})

test_that("read_annotations rejects empty files", {
  f <- tempfile()
  file.create(f)
  expect_error(read_annotations(f), "empty")
  unlink(f)
})

test_that("read_wfdb_record and prepare_record integrate with annotations", {
  d <- tempfile()
  dir.create(file.path(d, "atr"), recursive = TRUE)
  write_header(d, "001", 16L, n_sig = 1L, n_samples = 200L, desc = "MLII")
  con <- file(file.path(d, "001.dat"), "wb")
  sig <- as.integer(round(100 * sin(seq_len(200) / 5)))
  writeBin(word_bytes(sig + 500L), con)
  close(con)
  awords <- c(1L * 1024L + 50L, 1L * 1024L + 50L, 5L * 1024L + 50L,
              1L * 1024L + 40L)
  con <- file(file.path(d, "atr", "001.atr"), "wb")
  writeBin(word_bytes(awords), con)
  close(con)
  ann <- read_annotations(file.path(d, "atr", "001.atr"))
  expect_equal(ann$sample, c(50L, 100L, 150L, 190L))
  raw_sig <- read_wfdb("001", dir = d)
  X <- prepare_record(raw_sig$MLII, ann, record = "001", wsize = 5,
                      drop_first_last = FALSE)
  expect_equal(nrow(X), 4L)
  expect_equal(attr(X, "labels"), c("N", "N", "V", "N"))
  unlink(d, recursive = TRUE)
})

test_that("read_database reads every complete record and keeps the patient", {
  d <- tempfile()
  dir.create(d)
  for (r in c("001", "002")) {
    write_header(d, r, 16L, n_sig = 1L, n_samples = 200L, desc = "II")
    if (r == "002") cat("# <age>: 50 <sex>: F patient 7\n", file = file.path(d, "002.hea"),
                        append = TRUE)
    con <- file(file.path(d, paste0(r, ".dat")), "wb")
    writeBin(word_bytes(as.integer(round(100 * sin(seq_len(200) / 5))) + 500L), con)
    close(con)
    con <- file(file.path(d, paste0(r, ".atr")), "wb")
    writeBin(word_bytes(c(1L * 1024L + 50L, 1L * 1024L + 50L, 5L * 1024L + 50L,
                          1L * 1024L + 40L)), con)
    close(con)
  }
  write_header(d, "003", 16L, n_sig = 1L, n_samples = 200L, desc = "II")
  X <- read_database(d, wsize = 5L, lead = "II", min_beats = 1L)
  expect_named(X, c("001", "002"))
  expect_equal(ncol(X[[1]]), 11L)
  expect_null(attr(X[["001"]], "patient"))
  expect_equal(attr(X[["002"]], "patient"), "7")
  expect_length(read_database(d, wsize = 5L, min_beats = 10L), 0L)
  expect_error(read_database(d, wsize = 5L, lead = "V5", min_beats = 1L), "lead not found")
  unlink(d, recursive = TRUE)
})
