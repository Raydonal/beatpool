#' Colour palette for pooling methods
#'
#' A colour-blind friendly palette (Okabe-Ito) used consistently in all
#' figures.
#'
#' @return A named character vector of colours.
#' @export
method_palette <- function() {
  c("PCA" = "#000000", "Global PCA" = "#7F7F7F", "Tbar" = "#E69F00",
    "That" = "#56B4E9", "Procrustes" = "#009E73", "GPA" = "#0072B2",
    "DHmT" = "#CC79A7", "DCT" = "#D55E00", "Naive" = "#999999")
}

#' Consistent ggplot2 theme
#'
#' @param base_size Base font size.
#' @return A ggplot2 theme.
#' @export
theme_bp <- function(base_size = 12) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting", call. = FALSE)
  }
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
      axis.line = ggplot2::element_line(colour = "grey30", linewidth = 0.3),
      axis.ticks = ggplot2::element_line(colour = "grey30", linewidth = 0.3),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2),
      plot.subtitle = ggplot2::element_text(colour = "grey35")
    )
}

#' Plot communication cost against reconstruction error
#'
#' Places a set of operators in the cost--error plane. The cost is read from
#' the \code{numbers} attribute of a \code{bp_poolall} object, or supplied
#' directly, and the error is a named vector indexed by operator.
#'
#' @param x A \code{bp_poolall} object or a named numeric vector of costs.
#' @param error A named numeric vector of reconstruction errors.
#' @param labels Optional labels, one per operator. Defaults to the names.
#' @param k Number of coefficients, used in the axis label only.
#' @return A ggplot object.
#' @examples
#' sim <- simulate_beats(H = 6, n = 120, p = 20, q = 3, seed = 2)
#' all <- pool_all(sim$X, q = 3)
#' err <- vapply(all, function(T) mean(vapply(sim$X, reconstruction_error,
#'   numeric(1), T = T, k = 3)), numeric(1))
#' plot_cost_error(all, err)
#' @export
plot_cost_error <- function(x, error, labels = names(error), k = 3L) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting", call. = FALSE)
  }
  cost <- if (inherits(x, "bp_poolall")) attr(x, "numbers") else x
  if (is.null(names(error))) stop("error must be a named vector", call. = FALSE)
  common <- intersect(names(cost), names(error))
  if (!length(common)) {
    stop("no operator names in common between cost and error", call. = FALSE)
  }
  d <- data.frame(operator = common, numbers = as.numeric(cost[common]),
                  error = as.numeric(error[common]))
  d$label <- labels[match(common, names(error))]
  d$label[is.na(d$label)] <- common[is.na(d$label)]
  ggplot2::ggplot(d, ggplot2::aes(x = .data$numbers, y = .data$error)) +
    ggplot2::geom_point(colour = "#D55E00", size = 2.4) +
    ggplot2::geom_text(ggplot2::aes(label = .data$label), size = 2.9,
                       vjust = -0.9, check_overlap = TRUE) +
    ggplot2::scale_x_continuous(trans = "log1p",
                                breaks = c(0, 200, 500, 1000, 2000)) +
    ggplot2::labs(x = "Values transmitted per record (log scale)",
                  y = paste0("Held-out relative error at k = ", k)) +
    theme_bp()
}

#' Plot a sample of heartbeats
#'
#' @param X A beat matrix.
#' @param n Maximum number of beats to draw.
#' @param colour Optional vector of labels used to colour the beats.
#' @param alpha Line transparency.
#' @return A ggplot object.
#' @export
plot_beats <- function(X, n = 100L, colour = NULL, alpha = 0.25) {
  X <- as.matrix(X)
  idx <- if (nrow(X) > n) sample.int(nrow(X), n) else seq_len(nrow(X))
  df <- data.frame(
    time = rep(seq_len(ncol(X)), times = length(idx)),
    value = as.vector(t(X[idx, , drop = FALSE])),
    beat = rep(idx, each = ncol(X))
  )
  if (!is.null(colour)) df$group <- rep(colour[idx], each = ncol(X))
  p <- ggplot2::ggplot(df, ggplot2::aes(time, value))
  if (!is.null(colour)) {
    p <- p + ggplot2::geom_line(ggplot2::aes(group = beat, colour = group),
                                alpha = alpha)
  } else {
    p <- p + ggplot2::geom_line(ggplot2::aes(group = beat), alpha = alpha,
                                colour = "grey20")
  }
  p + ggplot2::labs(x = "Sample", y = "Standardised amplitude") + theme_bp()
}

#' Plot the mean heartbeat
#'
#' @param X A beat matrix.
#' @param ribbon Draw a pointwise confidence band.
#' @return A ggplot object.
#' @export
plot_mean_beat <- function(X, ribbon = TRUE) {
  X <- as.matrix(X)
  df <- data.frame(time = seq_len(ncol(X)), mean = colMeans(X),
                   sd = apply(X, 2L, stats::sd))
  p <- ggplot2::ggplot(df, ggplot2::aes(time, mean))
  if (ribbon) {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = mean - sd,
                                               ymax = mean + sd),
                                  fill = "#009E73", alpha = 0.2)
  }
  p + ggplot2::geom_line(colour = "#009E73", linewidth = 0.8) +
    ggplot2::labs(x = "Sample", y = "Mean standardised amplitude") +
    theme_bp()
}

#' Plot basis coefficients
#'
#' @param T A basis matrix.
#' @param k Number of components to plot.
#' @return A ggplot object.
#' @export
plot_basis <- function(T, k = min(3L, ncol(T))) {
  df <- basis_frame(as.matrix(T)[, seq_len(k), drop = FALSE])
  df$component <- factor(df$component, levels = seq_len(k))
  ggplot2::ggplot(df, ggplot2::aes(time, value, colour = component)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::facet_wrap(~component, ncol = 1, scales = "free_y") +
    ggplot2::labs(x = "Time instant", y = "Coefficient", colour = "Component") +
    theme_bp() + ggplot2::theme(legend.position = "none")
}

#' Plot a scree or cumulative variance curve
#'
#' @param values A numeric vector or a list of numeric vectors.
#' @param cumulative Plot cumulative proportions.
#' @param labels Optional names for the list entries.
#' @return A ggplot object.
#' @export
plot_scree <- function(values, cumulative = TRUE, labels = NULL) {
  if (!is.list(values)) values <- list(values)
  if (is.null(labels)) {
    labels <- base::names(values)
    if (is.null(labels)) labels <- paste0("Basis ", seq_along(values))
  }
  df <- do.call(rbind, lapply(seq_along(values), function(i) {
    v <- values[[i]]
    y <- if (cumulative) cumsum(v) / sum(v) else v / sum(v)
    data.frame(component = seq_along(v), value = y, basis = labels[i])
  }))
  ggplot2::ggplot(df, ggplot2::aes(component, value, colour = basis)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::labs(x = "Component", y = if (cumulative)
      "Cumulative variance proportion" else "Variance proportion",
      colour = NULL) + theme_bp()
}

#' Plot reconstruction error curves
#'
#' @param df A data frame with columns \code{k}, \code{error} and
#'   \code{method}.
#' @param facet Optional faceting variable.
#' @return A ggplot object.
#' @export
plot_reconstruction <- function(df, facet = NULL) {
  p <- ggplot2::ggplot(df, ggplot2::aes(k, error, colour = method)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::scale_colour_manual(values = method_palette(),
                                 drop = FALSE) +
    ggplot2::labs(x = "Number of components (k)",
                  y = "Relative reconstruction error", colour = NULL) +
    theme_bp()
  if (!is.null(facet)) p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet)))
  p
}

#' Plot clustering accuracy curves
#'
#' @param df A data frame with columns \code{k}, \code{value},
#'   \code{method} and optionally \code{K}.
#' @param metric Label for the vertical axis.
#' @return A ggplot object.
#' @export
plot_clustering <- function(df, metric = "Adjusted Rand index") {
  p <- ggplot2::ggplot(df, ggplot2::aes(k, value, colour = method)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::scale_colour_manual(values = method_palette(), drop = FALSE) +
    ggplot2::labs(x = "Number of components (k)", y = metric, colour = NULL) +
    theme_bp()
  if ("K" %in% names(df)) p <- p + ggplot2::facet_wrap(~K)
  p
}

#' Plot jackknife stability
#'
#' @param df A data frame returned by [jackknife_stability()].
#' @param type Which distance to plot: \code{"relative"} or \code{"subspace"}.
#' @return A ggplot object.
#' @export
plot_stability <- function(df, type = c("relative", "subspace")) {
  type <- match.arg(type)
  ggplot2::ggplot(df, ggplot2::aes(record, .data[[type]])) +
    ggplot2::geom_hline(yintercept = mean(df[[type]]), linetype = 2,
                        colour = "#D55E00") +
    ggplot2::geom_line(colour = "grey40", linewidth = 0.5) +
    ggplot2::geom_point(colour = "#0072B2", size = 1.8) +
    ggplot2::labs(x = "Omitted record", y = if (type == "relative")
      "Relative Frobenius change" else "Subspace distance") + theme_bp()
}

#' Plot simulation results
#'
#' @param df A data frame.
#' @param x Column mapped to the horizontal axis.
#' @param y Column mapped to the vertical axis.
#' @param colour Column mapped to colour.
#' @param facet Optional faceting variable.
#' @return A ggplot object.
#' @export
plot_simulation <- function(df, x, y, colour = "method", facet = NULL) {
  p <- ggplot2::ggplot(df, ggplot2::aes(.data[[x]], .data[[y]],
                                        colour = .data[[colour]])) +
    ggplot2::stat_summary(fun = mean, geom = "line", linewidth = 0.7) +
    ggplot2::stat_summary(fun = mean, geom = "point", size = 1.4) +
    ggplot2::scale_colour_manual(values = method_palette(), drop = FALSE) +
    ggplot2::labs(x = x, y = y, colour = NULL) + theme_bp()
  if (!is.null(facet)) p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet)))
  p
}
