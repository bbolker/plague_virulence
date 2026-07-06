## Plot two-strain burnout scan results.
## Put this file in the same folder as the output folder, then run:
## source("plot_two_strain_burnout_v9.R")

get_script_dir <- function() {
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    if (!is.null(frames[[i]]$ofile)) {
      return(dirname(normalizePath(frames[[i]]$ofile, winslash = "/", mustWork = TRUE)))
    }
  }
  getwd()
}

script_dir <- get_script_dir()
output_dir <- file.path(script_dir, "output")
csv_file <- file.path(output_dir, "two_strain_burnout_scan.csv")
if (!file.exists(csv_file)) {
  stop("Could not find output/two_strain_burnout_scan.csv.", call. = FALSE)
}

scan_results <- read.csv(csv_file, stringsAsFactors = FALSE)
R01_values <- sort(unique(scan_results$R01))
R02_values <- sort(unique(scan_results$R02))

make_matrix <- function(value_name) {
  z <- matrix(NA_real_, nrow = length(R01_values), ncol = length(R02_values))
  r1_index <- match(round(scan_results$R01, 10), round(R01_values, 10))
  r2_index <- match(round(scan_results$R02, 10), round(R02_values, 10))
  z[cbind(r1_index, r2_index)] <- scan_results[[value_name]]
  z
}

format_prob_label <- function(x) {
  out <- character(length(x))
  pow10 <- x > 0 & abs(log10(x) - round(log10(x))) < 1e-10 & x <= 1e-1
  out[pow10] <- paste0("1e", round(log10(x[pow10])))
  regular <- out == ""
  out[regular] <- format(x[regular], trim = TRUE, digits = 2)
  out
}

make_split_scale <- function(z, cutoff = 0.1, prob_floor = 1e-12) {
  z_pos <- z[is.finite(z) & z > 0]
  if (length(z_pos) == 0L) stop("No positive finite values to plot.", call. = FALSE)
  floor_value <- 10^floor(log10(min(prob_floor, min(z_pos))))
  floor_value <- max(floor_value, .Machine$double.xmin)

  map_prob <- function(p) {
    p <- pmax(p, floor_value)
    out <- numeric(length(p))
    low <- p <= cutoff
    out[low] <- 0.5 * (log10(p[low]) - log10(floor_value)) / (log10(cutoff) - log10(floor_value))
    out[!low] <- 0.5 + 0.5 * (p[!low] - cutoff) / (1 - cutoff)
    pmin(pmax(out, 0), 1)
  }

  map_matrix <- function(p) {
    out <- matrix(NA_real_, nrow = nrow(p), ncol = ncol(p))
    ok <- is.finite(p) & p > 0
    if (any(ok)) out[ok] <- map_prob(p[ok])
    out
  }

  list(map_prob = map_prob, map_matrix = map_matrix, floor = floor_value, cutoff = cutoff)
}

plot_heatmap <- function(value_name, main_title, file_name) {
  z <- make_matrix(value_name)
  z_range <- range(z, finite = TRUE, na.rm = TRUE)
  if (!all(is.finite(z_range))) stop("No finite values to plot.", call. = FALSE)

  split_scale <- make_split_scale(z, cutoff = 0.1, prob_floor = 1e-12)
  z_plot <- split_scale$map_matrix(z)
  z_plot_range <- c(0, 1)

  cols <- colorRampPalette(c("#313695", "#4575b4", "#74add1", "#abd9e9", "#e0f3f8", "#ffffbf", "#fee090", "#fdae61", "#f46d43", "#d73027", "#a50026"))(180)
  breaks <- seq(z_plot_range[1], z_plot_range[2], length.out = length(cols) + 1L)

  contour_levels <- c(1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1e-1, 0.2, 0.4, 0.8, 1)
  contour_levels <- contour_levels[contour_levels >= z_range[1] & contour_levels <= z_range[2]]

  tick_probs <- c(split_scale$floor, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1e-1, 0.2, 0.4, 0.8, 1)
  tick_probs <- sort(unique(tick_probs))
  tick_probs <- tick_probs[tick_probs >= split_scale$floor & tick_probs <= 1]
  tick_at <- split_scale$map_prob(tick_probs)

  pdf_file <- file.path(output_dir, file_name)
  pdf(pdf_file, width = 13, height = 9.5)
  layout(matrix(c(1, 2), nrow = 1), widths = c(5.4, 0.9))

  par(mar = c(4.8, 4.8, 3.2, 1.0), pty = "s")
  image(
    x = R01_values,
    y = R02_values,
    z = z_plot,
    xlab = "R01",
    ylab = "R02",
    main = main_title,
    col = cols,
    breaks = breaks,
    asp = 1,
    useRaster = TRUE
  )
  if (length(contour_levels) > 0L && z_range[1] < z_range[2]) {
    contour(
      R01_values, R02_values, z,
      levels = contour_levels,
      labels = format_prob_label(contour_levels),
      add = TRUE,
      drawlabels = TRUE,
      col = "black",
      lwd = 1.8,
      labcex = 0.8
    )
  }

  par(mar = c(4.8, 0.5, 3.2, 4.2), pty = "m")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = z_plot_range)
  for (i in seq_along(cols)) {
    rect(0, breaks[i], 1, breaks[i + 1L], col = cols[i], border = NA)
  }
  if (length(tick_probs) > 0L) {
    axis(4, at = tick_at, labels = format_prob_label(tick_probs), las = 1)
  }
  abline(h = 0.5, lwd = 1)
  box()
  mtext("probability", side = 4, line = 3.3)

  dev.off()
  layout(1)
  invisible(pdf_file)
}

plot_file <- plot_heatmap(
  "surv1",
  "Strain 1 persistence probability: 1 - Q1",
  "heatmap_strain1_persistence.pdf"
)

cat("Saved plot to ", normalizePath(plot_file, winslash = "/", mustWork = FALSE), "\n", sep = "")
