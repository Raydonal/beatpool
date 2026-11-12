#' Download a PhysioNet database
#'
#' Downloads the header, signal and annotation files of every record listed in
#' the \code{RECORDS} file of a PhysioNet database. Files already present in
#' \code{dir} are kept.
#'
#' @param db Database slug, such as \code{"mitdb"}, \code{"svdb"} or
#'   \code{"incartdb"}.
#' @param dir Destination directory.
#' @param version Database version.
#' @param ext File extensions downloaded for each record.
#' @param records Optional subset of records.
#' @return The record names, invisibly.
#' @examples
#' \dontrun{
#' download_physionet("mitdb", "data/mitdb")
#' }
#' @export
download_physionet <- function(db, dir, version = "1.0.0",
                               ext = c("hea", "dat", "atr"), records = NULL) {
  base <- sprintf("https://physionet.org/files/%s/%s/", db, version)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  fetch <- function(url, f, tries = 5L) {
    for (i in seq_len(tries)) {
      ok <- tryCatch({
        utils::download.file(url, f, mode = "wb", quiet = TRUE)
        TRUE
      }, error = function(e) FALSE, warning = function(w) FALSE)
      if (ok) return(invisible(TRUE))
      unlink(f)
      Sys.sleep(2 * i)
    }
    stop("could not download ", url, call. = FALSE)
  }
  if (is.null(records)) {
    rf <- file.path(dir, "RECORDS")
    if (!file.exists(rf)) fetch(paste0(base, "RECORDS"), rf)
    records <- trimws(readLines(rf, warn = FALSE))
    records <- records[nzchar(records)]
  }
  for (r in records) for (e in ext) {
    f <- file.path(dir, paste0(r, ".", e))
    if (!file.exists(f)) fetch(paste0(base, r, ".", e), f)
  }
  invisible(records)
}

#' Beat matrices of a WFDB database
#'
#' Reads every record of a database stored in WFDB format, segments one
#' channel around the annotated beats and keeps records with at least
#' \code{min_beats} beats. These are the settings of the accompanying article:
#' the first signal with \code{wsize = 30} for the MIT-BIH Arrhythmia Database
#' and \code{wsize = 15} for the MIT-BIH Supraventricular Arrhythmia Database,
#' and lead II with \code{wsize = 21} for the St Petersburg INCART database.
#'
#' @param dir Directory with the \code{.hea}, \code{.dat} and \code{.atr} files.
#' @param wsize Half window in samples; a beat has \code{2 * wsize + 1} samples.
#' @param lead Signal description to read, or \code{NULL} for the first signal.
#' @param min_beats Minimum number of beats of a record.
#' @param normalize Passed to [prepare_record()].
#' @return A named list of beat matrices. For records whose header names a
#'   patient, the attribute \code{patient} holds it.
#' @examples
#' \dontrun{
#' download_physionet("svdb", "data/svdb")
#' X <- read_database("data/svdb", wsize = 15)
#' }
#' @export
read_database <- function(dir, wsize, lead = NULL, min_beats = 100L,
                          normalize = "segment") {
  records <- sort(sub("\\.hea$", "", list.files(dir, pattern = "\\.hea$")))
  complete <- vapply(records, function(r)
    all(file.exists(file.path(dir, paste0(r, c(".hea", ".dat", ".atr"))))),
    logical(1))
  records <- records[complete]
  read_one <- function(r) {
    hea <- file.path(dir, paste0(r, ".hea"))
    hdr <- read_wfdb_header(hea)
    desc <- trimws(hdr$signals$description)
    ch <- if (is.null(lead)) 1L else which(desc == lead)
    if (length(ch) != 1L) stop("lead not found in record ", r, call. = FALSE)
    sig <- read_wfdb(r, dir = dir, channels = ch)
    ann <- read_annotations(file.path(dir, paste0(r, ".atr")))
    X <- prepare_record(sig[[2L]], ann, record = r, fs = hdr$fs,
                        wsize = wsize, lead = desc[ch], normalize = normalize,
                        drop_first_last = TRUE)
    txt <- grep("patient", readLines(hea, warn = FALSE), value = TRUE)
    if (length(txt))
      attr(X, "patient") <- sub(".*patient\\s+([0-9]+).*", "\\1", txt[1])
    X
  }
  out <- lapply(records, read_one)
  names(out) <- records
  out[vapply(out, nrow, integer(1)) >= min_beats]
}
