## Plot the focused stochastic-validation figures.
## Run after smoke or validation scripts:
##   Rscript fadeout/logistic_burnout/stochastic_validation/plot_stochastic_validation.R

source(file.path(
  "fadeout", "logistic_burnout", "stochastic_validation",
  "stochastic_validation_functions.R"
))
.require_pkg("ggplot2")

module_dir <- file.path("fadeout", "logistic_burnout", "stochastic_validation")
opt <- parse_cli_options(list(
  use_smoke = FALSE,
  output_dir = file.path(module_dir, "outputs"),
  figure_dir = file.path(module_dir, "figures")
))
output_dir <- ensure_dir(as.character(opt$output_dir))
figure_dir <- ensure_dir(as.character(opt$figure_dir))

full_file <- file.path(output_dir, "full_stochastic_validation.csv")
engine_file <- file.path(output_dir, "simulation_engine_comparison.csv")
if (!file.exists(full_file) || as_logical(opt$use_smoke)) {
  full_file <- file.path(output_dir, "smoke_test_results.csv")
}
if (!file.exists(full_file)) stop("No validation table found")

full <- read.csv(full_file)
engine <- if (file.exists(engine_file)) read.csv(engine_file) else NULL
full$R0_minus_1 <- full$R0 - 1
full$K_factor <- factor(full$K)

save_plot <- function(plot, filename, width = 9, height = 7) {
  ggplot2::ggsave(file.path(figure_dir, paste0(filename, ".pdf")),
                  plot, width = width, height = height)
  ggplot2::ggsave(file.path(figure_dir, paste0(filename, ".png")),
                  plot, width = width, height = height, dpi = 180)
}

theme_validation <- function() {
  ggplot2::theme_bw() +
    ggplot2::theme(
      plot.subtitle = ggplot2::element_text(size = 9),
      plot.caption = ggplot2::element_text(hjust = 0, size = 7),
      plot.margin = ggplot2::margin(8, 8, 12, 8),
      legend.position = "right"
    )
}

main_r <- full[abs(full$r - 0.1) < 1e-12, ]
if (!nrow(main_r)) stop("No r = 0.1 rows found for the main validation figure")

main_plot <- ggplot2::ggplot(
  main_r, ggplot2::aes(R0_minus_1, P_uncond_sim)
) +
  ggplot2::geom_line(ggplot2::aes(
    y = P_uncond_approx,
    colour = "Analytical approximation"
  ), linewidth = 0.7) +
  ggplot2::geom_line(ggplot2::aes(
    y = p_est_approx,
    colour = "Early-establishment reference"
  ), linewidth = 0.4) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = P_uncond_sim_ci_low,
                 ymax = P_uncond_sim_ci_high),
    width = 0
  ) +
  ggplot2::geom_point(
    ggplot2::aes(colour = "Stochastic simulation"),
    size = 1.8
  ) +
  ggplot2::scale_colour_manual(
    name = NULL,
    values = c(
      "Stochastic simulation" = "#2c7fb8",
      "Analytical approximation" = "black",
      "Early-establishment reference" = "grey45"
    ),
    breaks = c(
      "Stochastic simulation",
      "Analytical approximation",
      "Early-establishment reference"
    )
  ) +
  ggplot2::facet_wrap(~K_factor) +
  ggplot2::scale_x_log10() +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(
    x = "R0 - 1 (log scale)",
    y = "Unconditional persistence probability",
    title = "Full stochastic validation of unconditional persistence",
    subtitle = "Non-seasonal single-patch logistic S-I model; r = 0.1; I0 = 1",
    caption = paste(
      "Points and intervals are stochastic adaptive tau-leap estimates",
      "with 95% Wilson Monte Carlo CI."
    )
  ) +
  theme_validation()
save_plot(main_plot, "validation_unconditional_probability_scale_legend")

cond_plot <- ggplot2::ggplot(
  main_r, ggplot2::aes(R0_minus_1, P_cond_sim_established)
) +
  ggplot2::geom_line(ggplot2::aes(
    y = P_cond_approx,
    colour = "Analytical approximation (round m)"
  ), linewidth = 0.7) +
  ggplot2::geom_line(ggplot2::aes(
    y = P_cond_approx_continuous,
    colour = "Analytical approximation (continuous m)"
  ), linewidth = 0.5, linetype = "22") +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = P_cond_sim_established_ci_low,
                 ymax = P_cond_sim_established_ci_high),
    width = 0
  ) +
  ggplot2::geom_point(
    ggplot2::aes(colour = "Stochastic simulation"),
    size = 1.8
  ) +
  ggplot2::scale_colour_manual(
    name = NULL,
    values = c(
      "Stochastic simulation" = "#2c7fb8",
      "Analytical approximation (round m)" = "black",
      "Analytical approximation (continuous m)" = "#d95f02"
    ),
    breaks = c(
      "Stochastic simulation",
      "Analytical approximation (round m)",
      "Analytical approximation (continuous m)"
    )
  ) +
  ggplot2::facet_wrap(~K_factor) +
  ggplot2::scale_x_log10() +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(
    x = "R0 - 1 (log scale)",
    y = "Persistence probability, conditional on not fizzling",
    title = "Stochastic validation conditional on escaping early fizzle",
    subtitle = "Non-seasonal single-patch logistic S-I model; r = 0.1; I0 = 1",
    caption = paste(
      "P_cond_sim_established = n_persistence / (n_total - n_fizzle).",
      "Points and intervals are stochastic adaptive tau-leap estimates",
      "with 95% Wilson Monte Carlo CI. The dashed orange curve uses the",
      "unrounded m = K*y* directly as the exponent (Eq. 22/26/27 of",
      "Parsons et al. 2024), instead of rounding m to the nearest integer",
      "(solid black curve, the default used elsewhere in this repository)."
    )
  ) +
  theme_validation()
save_plot(cond_plot, "validation_conditional_probability_round_vs_continuous")

if (!is.null(engine) && nrow(engine)) {
  engine$label <- ifelse(
    engine$engine == "tau_leap",
    paste0("fixed tau dt=", engine$dt),
    ifelse(engine$engine == "adaptive_tau", "adaptive tau", "Gillespie")
  )
  engine_plot <- ggplot2::ggplot(
    engine,
    ggplot2::aes(label, P_uncond_sim, colour = label)
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = P_uncond_sim_ci_low,
                   ymax = P_uncond_sim_ci_high),
      width = 0
    ) +
    ggplot2::geom_point(size = 2.3) +
    ggplot2::facet_wrap(~paste0("R0=", R0, ", K=", K), scales = "free_x") +
    ggplot2::labs(
      x = NULL,
      y = "Unconditional persistence probability",
      colour = NULL,
      title = "Simulator and timestep sensitivity",
      subtitle = "Non-seasonal single-patch logistic S-I model; representative cells; I0 = 1",
      caption = "Intervals are 95% Wilson Monte Carlo CI."
    ) +
    theme_validation() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  save_plot(engine_plot, "validation_simulator_sensitivity", 10, 7)
} else {
  message("Skipping simulator sensitivity plot: no engine comparison CSV found")
}

cat("Focused plots written to ", figure_dir, "\n", sep = "")
