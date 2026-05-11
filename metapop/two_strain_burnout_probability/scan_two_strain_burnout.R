## two-strain burnout parameter scan.
## source("scan_two_strain_burnout.R")

if (!requireNamespace("burnout", quietly = TRUE)) {
  stop("Please install the burnout package: remotes::install_github('davidearn/burnout')", call. = FALSE)
}
if (!requireNamespace("gsl", quietly = TRUE)) {
  stop("Please install the gsl package: install.packages('gsl')", call. = FALSE)
}

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
source(file.path(script_dir, "two_strain_burnout.R"))

output_dir <- file.path(script_dir, "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

## Parameters to scan.
R01_values <- seq(1.2, 3.0, by = 0.05)
R02_values <- seq(1.2, 3.0, by = 0.05)

## Parameters held fixed.
epsilon <- 0.01
N <- 1e6
k1 <- 1
k2 <- 1
dt <- 0.01
t_max <- 200
n_x <- 4000L
q_method <- "q_approx"

rows <- vector("list", length(R01_values) * length(R02_values))
idx <- 1L

for (R01 in R01_values) {
  for (R02 in R02_values) {
    cat("Running R01 =", R01, ", R02 =", R02, "\n")

    ans <- tryCatch(
      two_strain_wave_burnout(
        R01 = R01,
        R02 = R02,
        epsilon = epsilon,
        N = N,
        k1 = k1,
        k2 = k2,
        dt = dt,
        t_max = t_max,
        n_x = n_x,
        q_method = q_method,
        store_trajectory = FALSE
      ),
      error = function(e) e
    )

    if (inherits(ans, "error")) {
      rows[[idx]] <- data.frame(
        R01 = R01, R02 = R02,
        boundary_event = "error",
        entry_case = NA_character_,
        Q1 = NA_real_, Q2 = NA_real_,
        surv1 = NA_real_, surv2 = NA_real_,
        both_burnout = NA_real_,
        strain1_only_burnout = NA_real_,
        strain2_only_burnout = NA_real_,
        neither_burnout = NA_real_,
        error_message = conditionMessage(ans)
      )
    } else if (!isTRUE(ans$boundary_event_found)) {
      rows[[idx]] <- data.frame(
        R01 = R01, R02 = R02,
        boundary_event = "no",
        entry_case = NA_character_,
        Q1 = NA_real_, Q2 = NA_real_,
        surv1 = NA_real_, surv2 = NA_real_,
        both_burnout = NA_real_,
        strain1_only_burnout = NA_real_,
        strain2_only_burnout = NA_real_,
        neither_burnout = NA_real_,
        error_message = NA_character_
      )
    } else {
      Q1 <- unname(ans$strain_burnout["Q1"])
      Q2 <- unname(ans$strain_burnout["Q2"])
      out <- ans$outcomes_by_strain

      rows[[idx]] <- data.frame(
        R01 = R01, R02 = R02,
        boundary_event = ans$boundary_event,
        entry_case = ans$entry_case,
        Q1 = Q1,
        Q2 = Q2,
        surv1 = 1 - Q1,
        surv2 = 1 - Q2,
        both_burnout = unname(out["both_burnout"]),
        strain1_only_burnout = unname(out["strain1_only_burnout"]),
        strain2_only_burnout = unname(out["strain2_only_burnout"]),
        neither_burnout = unname(out["neither_burnout"]),
        error_message = NA_character_
      )
    }

    idx <- idx + 1L
  }
}

scan_results <- do.call(rbind, rows)
csv_file <- file.path(output_dir, "two_strain_burnout_scan.csv")
write.csv(scan_results, csv_file, row.names = FALSE)

make_matrix <- function(value_name) {
  z <- matrix(NA_real_, nrow = length(R01_values), ncol = length(R02_values))
  for (i in seq_along(R01_values)) {
    for (j in seq_along(R02_values)) {
      hit <- scan_results$R01 == R01_values[i] & scan_results$R02 == R02_values[j]
      z[i, j] <- scan_results[[value_name]][hit][1]
    }
  }
  z
}

format_prob_label <- function(x) {
  out <- character(length(x))
  out[x == 0] <- "0"
  pos <- x > 0
  pow10 <- pos & abs(log10(x) - round(log10(x))) < 1e-10 & x <= 1e-1
  out[pow10] <- paste0("1e", round(log10(x[pow10])))
  regular <- out == ""
  out[regular] <- format(x[regular], trim = TRUE, digits = 2)
  out
}

plot_heatmap <- function(value_name, main_title, file_name) {
  z <- make_matrix(value_name)
  z_range <- range(z, finite = TRUE, na.rm = TRUE)
  if (!all(is.finite(z_range))) stop("No finite values to plot.", call. = FALSE)

  prob_floor <- 1e-12
  z_plot <- log10(pmax(z, prob_floor))
  z_plot_range <- range(z_plot, finite = TRUE, na.rm = TRUE)
  if (!all(is.finite(z_plot_range))) stop("No finite log-scale values to plot.", call. = FALSE)
  if (z_plot_range[1] == z_plot_range[2]) {
    z_plot_range <- z_plot_range + c(-0.5, 0.5)
  }

  cols <- hcl.colors(120, "Viridis", rev = FALSE)
  breaks <- seq(z_plot_range[1], z_plot_range[2], length.out = length(cols) + 1L)

  contour_levels <- sort(unique(c(0, 10^seq(-12, -1, by = 1), seq(0.1, 1, by = 0.1))))
  contour_levels <- contour_levels[contour_levels >= z_range[1] & contour_levels <= z_range[2]]

  pdf_file <- file.path(output_dir, file_name)
  pdf(pdf_file, width = 12.5, height = 9.5)
  layout(matrix(c(1, 2), nrow = 1), widths = c(5.2, 0.8))

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
    asp = 1
  )
  if (length(contour_levels) > 0L && z_range[1] < z_range[2]) {
    contour(
      R01_values, R02_values, z,
      levels = contour_levels,
      labels = format_prob_label(contour_levels),
      add = TRUE,
      drawlabels = TRUE,
      col = "black",
      lwd = 0.65,
      labcex = 0.6
    )
  }

  par(mar = c(4.8, 0.5, 3.2, 4.2), pty = "m")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = z_plot_range)
  for (i in seq_along(cols)) {
    rect(0, breaks[i], 1, breaks[i + 1L], col = cols[i], border = NA)
  }
  tick_probs <- sort(unique(c(prob_floor, 10^seq(-12, -1, by = 1), 0.2, 0.5, 1)))
  tick_probs <- tick_probs[tick_probs >= 10^z_plot_range[1] & tick_probs <= 10^z_plot_range[2]]
  axis(4, at = log10(tick_probs), labels = format_prob_label(tick_probs), las = 1)
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

cat("\nSaved results to ", normalizePath(csv_file, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Saved plot to ", normalizePath(plot_file, winslash = "/", mustWork = FALSE), "\n", sep = "")
