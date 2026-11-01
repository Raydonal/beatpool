#' WFDB annotation codes
#'
#' A table mapping the numeric annotation codes used by WFDB/MIT-BIH files to
#' their conventional symbols and a description.
#'
#' @format A data frame with columns \code{code}, \code{symbol},
#'   \code{description}, \code{beat} and \code{normal}.
#' @export
ann_codes <- data.frame(
  code = c(0:14, 16:21, 23:30, 36),
  symbol = c("", "N", "L", "R", "a", "V", "F", "J", "A", "S", "E", "j",
             "/", "Q", "~", "s", "T", "*", "D", "\"", "@", "B", "[", "]",
             "e", "n", "?", "?", "r", "f"),
  description = c(
    "not an annotation", "normal beat", "left bundle branch block beat",
    "right bundle branch block beat", "aberrated atrial premature beat",
    "premature ventricular contraction", "fusion of ventricular and normal beat",
    "nodal (junctional) premature beat", "atrial premature beat",
    "supraventricular premature beat", "ventricular escape beat",
    "nodal (junctional) escape beat", "paced beat", "unclassifiable beat",
    "signal quality change", "ST segment change", "T-wave change", "systole",
    "diastole", "comment annotation", "measurement annotation",
    "bundle branch block beat", "start of ventricular flutter/fibrillation",
    "end of ventricular flutter/fibrillation", "atrial escape beat",
    "supraventricular escape beat", "link to external data",
    "link to external data", "R-on-T premature ventricular contraction",
    "fusion of paced and normal beat"),
  beat = c(FALSE, rep(TRUE, 13), FALSE, rep(FALSE, 6), TRUE, FALSE, FALSE,
           TRUE, TRUE, FALSE, FALSE, TRUE, TRUE),
  normal = c(FALSE, TRUE, rep(FALSE, 28)),
  stringsAsFactors = FALSE
)

#' Read a WFDB header file
#'
#' Parses the \code{.hea} file of a WFDB record. Formats 212, 16 and 80 are
#' supported by [read_wfdb()].
#'
#' @param file Path to a \code{.hea} file.
#' @return A list with the record name, sampling frequency, number of samples,
#'   and a data frame describing each signal.
#' @export
read_wfdb_header <- function(file) {
  lines <- readLines(file, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  if (!length(lines)) stop("empty WFDB header: ", file, call. = FALSE)
  first <- strsplit(trimws(lines[1]), "[[:space:]]+")[[1]]
  record <- first[1]
  n_sig <- as.integer(first[2])
  fs <- as.numeric(first[3])
  n_samples <- if (length(first) >= 4 && grepl("^[0-9]+$", first[4])) {
    as.integer(first[4])
  } else NA_integer_

  sig_lines <- lines[seq_len(n_sig) + 1L]
  parse_sig <- function(line) {
    f <- strsplit(trimws(line), "[[:space:]]+")[[1]]
    if (length(f) < 8) stop("malformed signal line: ", line, call. = FALSE)
    gain_str <- f[3]
    m <- regmatches(gain_str, regexec(
      "^([-+0-9.eE]+)(\\(([-+0-9]+)\\))?((/)(.*))?$", gain_str))[[1]]
    gain <- as.numeric(m[2])
    baseline <- if (length(m) >= 4 && nzchar(m[4])) as.numeric(m[4]) else NA_real_
    units <- if (length(m) >= 7) m[7] else ""
    data.frame(
      file = f[1], format = as.integer(f[2]), gain = gain,
      baseline = baseline, units = units,
      adc_res = as.integer(f[4]), adc_zero = as.integer(f[5]),
      init_value = as.numeric(f[6]), checksum = as.numeric(f[7]),
      block_size = as.integer(f[8]),
      description = if (length(f) > 8) paste(f[-(1:8)], collapse = " ") else "",
      stringsAsFactors = FALSE
    )
  }
  signals <- do.call(rbind, lapply(sig_lines, parse_sig))
  signals$baseline[is.na(signals$baseline)] <- signals$adc_zero[is.na(signals$baseline)]
  list(record = record, n_sig = n_sig, fs = fs, n_samples = n_samples,
       signals = signals, path = normalizePath(file, mustWork = FALSE))
}

#' Read a WFDB annotation file
#'
#' Reads MIT-BIH style annotation files. Each annotation is one 16-bit word
#' containing a 10-bit time increment and a 6-bit code. Auxiliary strings
#' (code 63) and skip records (code 59) are handled.
#'
#' @param file Path to a \code{.atr} (or other annotation) file.
#' @return A data frame with one row per annotation and columns
#'   \code{sample} (0-based WFDB sample), \code{index} (1-based sample),
#'   \code{code}, \code{symbol}, \code{beat} and \code{normal}.
#' @export
read_annotations <- function(file) {
  size <- file.size(file)
  if (is.na(size) || size < 2L) stop("annotation file is empty: ", file, call. = FALSE)
  raw <- readBin(file, "integer", n = size %/% 2L, size = 2L,
                 endian = "little", signed = FALSE)
  n <- length(raw)
  i <- 1L
  t <- 0L
  samples <- integer(0)
  codes <- integer(0)
  while (i <= n) {
    w <- raw[i]
    i <- i + 1L
    code <- w %/% 1024L
    val <- w %% 1024L
    if (code == 63L) {                     # AUX string
      if (val > 0L) i <- i + ceiling(val / 2)
      next
    }
    if (code == 59L) {                     # SKIP: 32-bit absolute time
      if (i + 1L > n) break
      skip_to <- as.numeric(raw[i]) + as.numeric(raw[i + 1L]) * 65536
      i <- i + 2L
      ## 0xFFFFFFFF and similar sentinels are not valid sample numbers
      if (is.finite(skip_to) && skip_to <= 2147483647) t <- skip_to
      next
    }
    if (code %in% 60:62) {                 # NUM / SUB / CHN
      i <- i + 1L
      next
    }
    t <- t + val
    if (code == 0L) next                  # not an annotation
    samples <- c(samples, t)
    codes <- c(codes, code)
  }
  idx <- match(codes, ann_codes$code)
  out <- data.frame(
    sample = samples,
    index = samples + 1L,
    code = codes,
    symbol = ifelse(is.na(idx), "?", ann_codes$symbol[idx]),
    description = ifelse(is.na(idx), "unknown annotation",
                         ann_codes$description[idx]),
    beat = ifelse(is.na(idx), FALSE, ann_codes$beat[idx]),
    normal = ifelse(is.na(idx), FALSE, ann_codes$normal[idx]),
    stringsAsFactors = FALSE
  )
  attr(out, "file") <- file
  out
}

.decode_212 <- function(bytes, n_sig) {
  n3 <- length(bytes) %/% 3L
  b <- bytes[seq_len(n3 * 3L)]
  b1 <- as.integer(b[seq(1L, length(b), by = 3L)])
  b2 <- as.integer(b[seq(2L, length(b), by = 3L)])
  b3 <- as.integer(b[seq(3L, length(b), by = 3L)])
  s1 <- b1 + bitwAnd(b2, 15L) * 256L
  s2 <- bitwShiftR(b2, 4L) + b3 * 16L
  out <- integer(2L * n3)
  out[seq(1L, length(out), by = 2L)] <- s1
  out[seq(2L, length(out), by = 2L)] <- s2
  out[out > 2047L] <- out[out > 2047L] - 4096L
  out
}

#' Read a WFDB signal record
#'
#' Pure-R reader for WFDB signal files in formats 212, 16 and 80. No external
#' WFDB software is required.
#'
#' @param record Record name (e.g. \code{"100"}).
#' @param dir Directory containing the \code{.hea} and \code{.dat} files.
#' @param channels Integer vector of channels to return; \code{NULL} returns all.
#' @param from First sample (1-based) to return.
#' @param to Last sample (1-based) to return; \code{NULL} returns to the end.
#' @return A data frame with a \code{time} column (seconds) and one column per
#'   channel, named after the header descriptions. The parsed header is stored
#'   in the \code{"header"} attribute.
#' @export
read_wfdb <- function(record, dir = ".", channels = NULL, from = 1L, to = NULL) {
  hea <- file.path(dir, paste0(record, ".hea"))
  header <- read_wfdb_header(hea)
  if (is.null(channels)) channels <- seq_len(header$n_sig)
  if (any(channels < 1L | channels > header$n_sig)) {
    stop("invalid channel selection", call. = FALSE)
  }
  dat_files <- unique(header$signals$file)
  if (length(dat_files) != 1L) {
    stop("records with multiple signal files are not supported", call. = FALSE)
  }
  dat <- file.path(dir, dat_files)
  formats <- unique(header$signals$format)
  if (length(formats) != 1L || !formats %in% c(212L, 16L, 80L)) {
    stop("unsupported WFDB format: ", paste(formats, collapse = ", "),
         call. = FALSE)
  }
  size <- file.size(dat)
  if (formats == 212L) {
    samples <- .decode_212(readBin(dat, "raw", n = size), header$n_sig)
  } else if (formats == 16L) {
    samples <- readBin(dat, "integer", n = size %/% 2L, size = 2L,
                       endian = "little", signed = TRUE)
  } else {
    samples <- as.integer(readBin(dat, "raw", n = size))
  }
  n_total <- length(samples) %/% header$n_sig
  if (is.na(header$n_samples)) header$n_samples <- n_total
  if (is.null(to)) to <- n_total
  from <- max(1L, as.integer(from))
  to <- min(as.integer(to), n_total)
  if (from > to) stop("empty sample range", call. = FALSE)
  sel <- seq.int(from, to)
  out <- data.frame(time = (sel - 1L) / header$fs)
  for (ch in channels) {
    sig <- header$signals[ch, ]
    v <- samples[seq.int(ch, length(samples), by = header$n_sig)][sel]
    ## some records declare gain 0 when the units are unspecified
    gain <- if (is.finite(sig$gain) && sig$gain != 0) sig$gain else 1
    v <- (v - sig$baseline) / gain
    nm <- if (nzchar(sig$description)) sig$description else paste0("ch", ch)
    out[[nm]] <- v
  }
  attr(out, "header") <- header
  attr(out, "channels") <- channels
  out
}

#' Read a MIT-BIH CSV exported by \code{rdsamp}
#'
#' @param file Path to the \code{.csv} file.
#' @param lead One of \code{"MLII"} or \code{"V5"}.
#' @return A data frame with \code{time} and \code{x}.
#' @export
read_mitdb_csv <- function(file, lead = c("MLII", "V5")) {
  lead <- match.arg(lead)
  d <- readr::read_csv(file, col_names = FALSE, skip = 2,
                       show_col_types = FALSE, progress = FALSE)
  names(d) <- c("time", "MLII", "V5")
  data.frame(time = d$time, x = d[[lead]])
}
