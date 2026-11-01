#' Segment ECG heartbeats around R peaks
#'
#' Extracts fixed-length windows around annotated R peaks. Each window has
#' \code{2 * wsize + 1} samples. By default the window is standardised with the
#' mean and standard deviation of the whole beat segment delimited by the
#' midpoints between consecutive R peaks, which is the preprocessing used in the
#' original MIT-BIH analysis.
#'
#' @param x Numeric ECG signal.
#' @param peaks Integer vector of R-peak positions (1-based).
#' @param labels Optional beat labels, one per peak.
#' @param normal Optional logical vector, one per peak.
#' @param record Record identifier stored in the result.
#' @param fs Sampling frequency in Hz.
#' @param wsize Half-window size. The window length is \code{2 * wsize + 1}.
#' @param normalize One of \code{"segment"} (default), \code{"qrs"} or
#'   \code{"none"}.
#' @param drop_edges Drop beats whose window falls outside the signal.
#' @param drop_first_last Drop the first and last annotated beats.
#' @return A matrix of beats (rows) by time samples (columns), with attributes
#'   \code{peaks}, \code{labels}, \code{normal}, \code{record}, \code{fs},
#'   \code{wsize} and \code{normalize}.
#' @export
segment_beats <- function(x, peaks, labels = NULL, normal = NULL,
                          record = "record", fs = 360, wsize = 30L,
                          normalize = c("segment", "qrs", "none"),
                          drop_edges = TRUE, drop_first_last = FALSE) {
  normalize <- match.arg(normalize)
  peaks <- as.integer(peaks)
  if (is.null(labels)) labels <- rep(NA_character_, length(peaks))
  labels <- as.character(labels)
  if (is.null(normal)) normal <- rep(NA, length(peaks))
  normal <- as.logical(normal)
  if (length(labels) != length(peaks) || length(normal) != length(peaks)) {
    stop("peaks, labels and normal must have the same length", call. = FALSE)
  }
  ord <- order(peaks)
  peaks <- peaks[ord]; labels <- labels[ord]; normal <- normal[ord]
  keep_dup <- !duplicated(peaks)
  peaks <- peaks[keep_dup]; labels <- labels[keep_dup]; normal <- normal[keep_dup]

  n <- length(x)
  p <- 2L * as.integer(wsize) + 1L
  half <- as.integer(wsize)
  valid <- peaks - half >= 1L & peaks + half <= n
  if (drop_first_last && length(peaks) > 2L) {
    valid[c(1L, length(peaks))] <- FALSE
  }
  if (drop_edges) {
    peaks <- peaks[valid]; labels <- labels[valid]; normal <- normal[valid]
  } else if (any(!valid)) {
    warning(sum(!valid), " beat(s) outside the signal were dropped",
            call. = FALSE)
    peaks <- peaks[valid]; labels <- labels[valid]; normal <- normal[valid]
  }

  m <- length(peaks)
  if (!m) {
    out <- matrix(numeric(0), 0L, p)
    attr(out, "peaks") <- integer(0); attr(out, "labels") <- character(0)
    attr(out, "normal") <- logical(0); attr(out, "record") <- record
    attr(out, "fs") <- fs; attr(out, "wsize") <- as.integer(wsize)
    attr(out, "normalize") <- normalize
    return(out)
  }

  X <- matrix(0, m, p)
  for (i in seq_len(m)) {
    r <- peaks[i]
    win <- x[(r - half):(r + half)]
    if (normalize == "none") {
      X[i, ] <- win
      next
    }
    if (normalize == "qrs") {
      mu <- mean(win); s <- stats::sd(win)
    } else {
      left <- if (i == 1L) 0L else floor((peaks[i - 1L] + peaks[i]) / 2)
      right <- if (i == m) n else floor((peaks[i] + peaks[i + 1L]) / 2)
      if (right - left < 2L) { mu <- mean(win); s <- stats::sd(win) }
      else { seg <- x[(left + 1L):right]; mu <- mean(seg); s <- stats::sd(seg) }
    }
    if (!is.finite(s) || s < .Machine$double.eps) s <- 1
    X[i, ] <- (win - mu) / s
  }
  attr(X, "peaks") <- peaks
  attr(X, "labels") <- labels
  attr(X, "normal") <- normal
  attr(X, "record") <- record
  attr(X, "fs") <- fs
  attr(X, "wsize") <- as.integer(wsize)
  attr(X, "normalize") <- normalize
  X
}

#' Prepare a record from a signal and its annotations
#'
#' @param x Numeric ECG signal.
#' @param ann Annotation data frame produced by [read_annotations()].
#' @param record Record identifier.
#' @param fs Sampling frequency.
#' @param wsize Half-window size.
#' @param lead Lead name.
#' @param normalize Normalisation method passed to [segment_beats()].
#' @param drop_first_last Drop the first and last annotated beats.
#' @return A matrix of beats with attributes, as returned by
#'   [segment_beats()].
#' @export
prepare_record <- function(x, ann, record = "record", fs = 360, wsize = 30L,
                           lead = "MLII",
                           normalize = c("segment", "qrs", "none"),
                           drop_first_last = TRUE) {
  normalize <- match.arg(normalize)
  ann <- ann[ann$beat, , drop = FALSE]
  ann <- ann[order(ann$index), , drop = FALSE]
  X <- segment_beats(
    x = x, peaks = ann$index, labels = ann$symbol, normal = ann$normal,
    record = record, fs = fs, wsize = wsize, normalize = normalize,
    drop_edges = TRUE, drop_first_last = drop_first_last
  )
  attr(X, "lead") <- lead
  X
}

#' Read and prepare a MIT-BIH record from CSV and annotation files
#'
#' @param record Record name, e.g. \code{"100"}.
#' @param data_dir Directory containing \code{<record>.csv} and
#'   \code{atr/<record>.atr}.
#' @param lead ECG lead to use.
#' @param wsize Half-window size.
#' @param normalize Normalisation method.
#' @param drop_first_last Drop the first and last annotated beats.
#' @return A matrix of beats with attributes.
#' @export
read_mitdb_record <- function(record, data_dir = ".", lead = c("MLII", "V5"),
                              wsize = 30L, normalize = c("segment", "qrs", "none"),
                              drop_first_last = TRUE) {
  lead <- match.arg(lead); normalize <- match.arg(normalize)
  csv <- file.path(data_dir, paste0(record, ".csv"))
  atr <- file.path(data_dir, "atr", paste0(record, ".atr"))
  if (!file.exists(csv)) stop("signal file not found: ", csv, call. = FALSE)
  if (!file.exists(atr)) stop("annotation file not found: ", atr, call. = FALSE)
  sig <- read_mitdb_csv(csv, lead = lead)
  ann <- read_annotations(atr)
  prepare_record(sig$x, ann, record = record, fs = 360, wsize = wsize,
                 lead = lead, normalize = normalize,
                 drop_first_last = drop_first_last)
}

#' Beat-level summaries
#'
#' @param X A beat matrix.
#' @return A data frame with one row per beat and columns \code{record},
#'   \code{beat}, \code{label}, \code{normal}, \code{mean}, \code{sd},
#'   \code{min}, \code{max} and \code{energy}.
#' @export
beat_summary <- function(X) {
  if (!is.matrix(X)) stop("X must be a matrix", call. = FALSE)
  data.frame(
    record = attr(X, "record") %||% NA_character_,
    beat = seq_len(nrow(X)),
    label = attr(X, "labels") %||% rep(NA_character_, nrow(X)),
    normal = attr(X, "normal") %||% rep(NA, nrow(X)),
    mean = rowMeans(X), sd = apply(X, 1L, stats::sd),
    min = apply(X, 1L, min), max = apply(X, 1L, max),
    energy = rowSums(X^2), stringsAsFactors = FALSE
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a
