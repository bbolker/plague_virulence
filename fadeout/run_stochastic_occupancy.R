library(plagueMetapop)
library(dplyr)
library(ggplot2)
library(here)

source(here::here("fadeout", "occupancy_functions.R"))
theme_set(theme_bw())

## One-factor-at-a-time exploration. Initial conditions and occupancy analysis
## follow fadeout/single_strain_metapop.R. Full occupancy after the initial
## decline is treated as an absorbing boundary for the displayed trajectory.
baseline <- list(
  R0 = 2.5,
  K = 1e4,
  r = 0.125,
  alpha = 1e-4,
  gamma = 1,
  n_patch = 200,
  dt = 0.1,
  t_max = 2000,
  I_outbreak = 10,
  initialization = "virgin_soil",
  initialization_seed = 101,
  simulation_seed = 102
)

parameter_values <- list(
  R0 = c(1.5, 2, 2.5, 3, 4, 5),
  K = c(1000, 3000, 10000, 30000, 100000),
  alpha = c(1e-5, 3e-5, 1e-4, 3e-4, 1e-3)
)

make_scenarios <- function(baseline, parameter_values) {
  unlist(lapply(names(parameter_values), function(parameter) {
    lapply(parameter_values[[parameter]], function(value) {
      params <- baseline
      params[[parameter]] <- value
      list(
        varied_parameter = parameter,
        parameter_value = value,
        params = params
      )
    })
  }), recursive = FALSE)
}

format_parameter <- function(name, value) {
  if (name %in% c("alpha")) {
    sprintf("%s = %g", name, value)
  } else if (name == "K") {
    sprintf("%s = %s", name, format(value, scientific = FALSE, big.mark = ","))
  } else {
    sprintf("%s = %g", name, value)
  }
}

fixed_parameter_subtitle <- function(varied_parameter, baseline) {
  fixed <- setdiff(c("R0", "K", "r", "alpha"), varied_parameter)
  paste(
    c(
      vapply(fixed, function(x) format_parameter(x, baseline[[x]]), character(1)),
      sprintf("n_patch = %d", baseline$n_patch),
      sprintf("S(0) = K - %d; I(0) = %d per patch",
              baseline$I_outbreak, baseline$I_outbreak),
      sprintf("dt = %g", baseline$dt),
      sprintf("t_max = %g", baseline$t_max)
    ),
    collapse = "; "
  )
}

outdir <- here::here("fadeout", "output", "stochastic_patch_occupancy")
figdir <- file.path(outdir, "figures")
datadir <- file.path(outdir, "data")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(datadir, recursive = TRUE, showWarnings = FALSE)

scenarios <- make_scenarios(baseline, parameter_values)
plot_only <- "--plot-only" %in% commandArgs(trailingOnly = TRUE)
if (plot_only) {
  result_file <- file.path(datadir, "stochastic_patch_occupancy_results.rds")
  if (!file.exists(result_file)) stop("Existing occupancy results not found: ", result_file)
  results <- readRDS(result_file)
  results <- Filter(function(x) x$varied_parameter %in% names(parameter_values), results)
  message("Plot-only mode: reusing existing R0, K, and alpha trajectories")
} else {
  gen <- compile_odin("euler_odin_def.R")
  results <- vector("list", length(scenarios))

  for (i in seq_along(scenarios)) {
    scenario <- scenarios[[i]]
    message(sprintf(
      "Running %s = %g (%d/%d)",
      scenario$varied_parameter,
      scenario$parameter_value,
      i,
      length(scenarios)
    ))

    result <- run_fadeout_occupancy(scenario$params, gen = gen)
    displayed <- apply_full_occupancy_absorbing_boundary(result$meta_summary)
    absorption_step <- attr(displayed, "absorption_step")

    results[[i]] <- list(
      varied_parameter = scenario$varied_parameter,
      parameter_value = scenario$parameter_value,
      params = scenario$params,
      raw_meta_summary = result$meta_summary,
      displayed_meta_summary = displayed,
      absorption_step = absorption_step,
      raw_curve_summary = result$curve_summary
    )
  }
}

occupancy_data <- bind_rows(lapply(results, function(x) {
  x$displayed_meta_summary |>
    select(step, occupied_patches, occupancy_fraction) |>
    mutate(
      varied_parameter = x$varied_parameter,
      parameter_value = x$parameter_value,
      parameter_label = format_parameter(
        x$varied_parameter, x$parameter_value
      )
    )
}))

curve_summaries <- bind_rows(lapply(results, function(x) {
  x$raw_curve_summary |>
    mutate(
      varied_parameter = x$varied_parameter,
      parameter_value = x$parameter_value,
      R0 = x$params$R0,
      K = x$params$K,
      r = x$params$r,
      alpha = x$params$alpha,
      absorption_step = x$absorption_step,
      absorbed_by_tmax = !is.na(x$absorption_step),
      displayed_final_occupancy =
        tail(x$displayed_meta_summary$occupancy_fraction, 1),
      .before = 1
    )
}))

palette <- c(
  "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"
)

for (parameter in names(parameter_values)) {
  plot_data <- occupancy_data |>
    filter(varied_parameter == parameter) |>
    mutate(parameter_label = factor(
      parameter_label,
      levels = vapply(
        parameter_values[[parameter]],
        function(value) format_parameter(parameter, value),
        character(1)
      )
    ))

  p_compare <- ggplot(
    plot_data,
    aes(step, occupancy_fraction, colour = parameter_label)
  ) +
    geom_line(linewidth = 0.65) +
    scale_colour_manual(values = palette[seq_along(parameter_values[[parameter]])]) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      labels = scales::percent_format(accuracy = 1)
    ) +
    labs(
      x = "Time (disease generations)",
      y = "Fraction of occupied patches",
      colour = sprintf("Varied parameter: %s", parameter),
      title = sprintf("Metapopulation occupancy: varying %s", parameter),
      subtitle = fixed_parameter_subtitle(parameter, baseline),
      caption = "After recovery to 100%, occupancy is treated as absorbing at 100%."
    ) +
    theme(
      legend.position = "right",
      plot.subtitle = element_text(size = 9),
      plot.caption = element_text(hjust = 0)
    )

  ggsave(
    file.path(figdir, sprintf("stochastic_patch_occupancy_compare_%s.pdf", parameter)),
    p_compare,
    width = 10,
    height = 6
  )
}

write.csv(
  occupancy_data,
  file.path(datadir, "stochastic_patch_occupancy_curves.csv"),
  row.names = FALSE
)
write.csv(
  curve_summaries,
  file.path(datadir, "stochastic_patch_occupancy_summaries.csv"),
  row.names = FALSE
)
saveRDS(results, file.path(datadir, "stochastic_patch_occupancy_results.rds"))

message("Finished. Outputs saved in: ", outdir)
