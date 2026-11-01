#' beatpool: pooled PCA bases for ECG heartbeats
#'
#' The package provides tools to segment heartbeats, estimate
#' patient-specific principal component bases, align and pool them across
#' records, compress and cluster beats, assess stability, and simulate
#' heterogeneous heartbeat populations.
#'
#' @keywords internal
#' @importFrom rlang .data
"_PACKAGE"

utils::globalVariables(c("time", "value", "beat", "group", "mean", "sd",
                         "component", "basis", "k", "error", "method",
                         "record"))
