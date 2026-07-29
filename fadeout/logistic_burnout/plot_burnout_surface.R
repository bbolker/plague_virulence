## Generate the exploratory R0-r surface for post-epidemic logistic burnout.
## Run from the repository root:
##   Rscript fadeout/logistic_burnout/plot_burnout_surface.R

source(file.path(
  "fadeout", "logistic_burnout", "logistic_burnout_functions.R"
))
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required")
}

module_dir <- file.path("fadeout", "logistic_burnout")
output_dir <- file.path(module_dir, "outputs")
figure_dir <- file.path(module_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

## K=10,000 is the baseline carrying capacity used repeatedly in the existing
## fadeout occupancy comparisons. I0=1 represents introduction by one infected
## host, while the boundary-entry lineage count is K*y_BL and is usually >1.
K_default <- 10000
I0_default <- 1
R0_values <- seq(1.05, 5, length.out = 41)
r_values <- seq(0.01, 0.5, length.out = 41)
grid <- expand.grid(R0 = R0_values, r = r_values)

start_time <- proc.time()[["elapsed"]]
results <- lapply(seq_len(nrow(grid)), function(i) {
  if (i %% 100 == 0 || i == 1L || i == nrow(grid)) {
    message("Computing grid cell ", i, "/", nrow(grid))
  }
  p <- grid[i, ]
  ans <- tryCatch(
    logistic_burnout_probability(
      R0 = p$R0, r = p$r, K = K_default, I0 = I0_default,
      lineage_count_method = "round",
      initial_tmax = 50, maximum_tmax = 1600, dt = 0.05,
      rtol = 1e-8, atol = 1e-10,
      rel.tol = 1e-7, subdivisions = 500
    ),
    error = function(e) list(
      R0 = p$R0, r = p$r, K = K_default, I0 = I0_default,
      x_in = NA_real_, y_BL = NA_real_,
      I_in = NA_real_, q1 = NA_real_, m_raw = NA_real_,
      m_used = NA_real_, P_burnout = NA_real_,
      status = paste0("error: ", conditionMessage(e)),
      boundary_entry_found = FALSE, integration_converged = FALSE,
      integration_absolute_error = NA_real_,
      integration_subdivisions = NA_integer_,
      t_peak = NA_real_, I_peak = NA_real_, t_in = NA_real_
    )
  )
  data.frame(
    R0 = ans$R0, r = ans$r, K = ans$K, I0 = ans$I0,
    x_in = ans$x_in, y_BL = ans$y_BL,
    I_in = ans$I_in, q1 = ans$q1, m_raw = ans$m_raw,
    m_used = ans$m_used, P_burnout = ans$P_burnout,
    status = ans$status,
    boundary_entry_found = ans$boundary_entry_found,
    integration_converged = ans$integration_converged,
    integration_absolute_error = ans$integration_absolute_error,
    integration_subdivisions = ans$integration_subdivisions,
    t_peak = ans$t_peak, I_peak = ans$I_peak, t_in = ans$t_in,
    stringsAsFactors = FALSE
  )
})
surface <- do.call(rbind, results)
runtime_seconds <- proc.time()[["elapsed"]] - start_time
write.csv(
  surface,
  file.path(output_dir, "logistic_burnout_R0_r_grid.csv"),
  row.names = FALSE
)

successful <- surface$boundary_entry_found &
  surface$integration_converged &
  is.finite(surface$P_burnout)
plot_surface <- surface
plot_surface$P_burnout[!successful] <- NA_real_

subtitle <- paste0(
  "Post-epidemic boundary-layer probability; K = ",
  format(K_default, big.mark = ","), "; I(0) = ", I0_default,
  "; m = max[1, round(K y_BL)]"
)
caption <- paste0(
  "The deterministic major-epidemic path starts at S(0)=K-1, I(0)=1. ",
  "\nThis conditional post-epidemic probability excludes early stochastic ",
  "extinction. Grey cells denote failed calculations."
)

p_burnout <- ggplot2::ggplot(
  plot_surface, ggplot2::aes(R0, r, fill = P_burnout)
) +
  ggplot2::geom_tile() +
  ggplot2::geom_contour(
    data = plot_surface[successful, ],
    ggplot2::aes(x = R0, y = r, z = P_burnout),
    inherit.aes = FALSE,
    breaks = c(0.1, 0.25, 0.5, 0.75, 0.9),
    colour = "white", linewidth = 0.3, alpha = 0.8
  ) +
  ggplot2::scale_fill_viridis_c(
    limits = c(0, 1), option = "magma", direction = -1,
    na.value = "grey70", name = "Burnout probability"
  ) +
  ggplot2::scale_x_continuous(expand = c(0, 0)) +
  ggplot2::scale_y_continuous(expand = c(0, 0)) +
  ggplot2::labs(
    x = "R0", y = "Logistic growth rate r per disease generation",
    title = "Post-epidemic burnout under logistic susceptible recovery",
    subtitle = subtitle, caption = caption
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_text(size = 9),
    plot.caption = ggplot2::element_text(hjust = 0, size = 8)
  )
ggplot2::ggsave(
  file.path(figure_dir, "logistic_burnout_R0_r_heatmap.png"),
  p_burnout, width = 9, height = 7, dpi = 180
)
ggplot2::ggsave(
  file.path(figure_dir, "logistic_burnout_R0_r_heatmap.pdf"),
  p_burnout, width = 9, height = 7
)

p_xin <- ggplot2::ggplot(
  surface,
  ggplot2::aes(R0, r, fill = ifelse(successful, x_in, NA_real_))
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_viridis_c(
    limits = c(0, 1), option = "viridis", na.value = "grey70",
    name = "Boundary-entry x_in"
  ) +
  ggplot2::scale_x_continuous(expand = c(0, 0)) +
  ggplot2::scale_y_continuous(expand = c(0, 0)) +
  ggplot2::labs(
    x = "R0", y = "Logistic growth rate r per disease generation",
    title = "Susceptible fraction at boundary-layer entry",
    subtitle = paste0(
      "First post-peak downward crossing of y = y_BL; ",
      "y_BL = r(R0-1)/R0^2; K = ",
      format(K_default, big.mark = ",")
    ),
    caption = paste0(
      "Grey cells denote failure to find the interpolated boundary\n",
      "entry or failure of the lineage integral."
    )
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_text(size = 9),
    plot.caption = ggplot2::element_text(hjust = 0, size = 8)
  )
ggplot2::ggsave(
  file.path(figure_dir, "logistic_boundary_xin_heatmap.png"),
  p_xin, width = 9, height = 7, dpi = 180
)
ggplot2::ggsave(
  file.path(figure_dir, "logistic_boundary_xin_heatmap.pdf"),
  p_xin, width = 9, height = 7
)

failure_table <- as.data.frame(table(
  status = surface$status, useNA = "ifany"
))
write.csv(
  failure_table,
  file.path(output_dir, "logistic_burnout_status_summary.csv"),
  row.names = FALSE
)

cat("Grid cells: ", nrow(surface), "\n", sep = "")
cat("Successful calculations: ", sum(successful), "/", nrow(surface),
    "\n", sep = "")
cat("Runtime (seconds): ", round(runtime_seconds, 1), "\n", sep = "")
cat("Burnout probability range: ",
    paste(signif(range(surface$P_burnout[successful]), 5), collapse = " to "),
    "\n", sep = "")
print(failure_table)
