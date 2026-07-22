## Generate tau=50 established occupancy and episode-classified transient,
## persistent, and censored occupancy/source-pressure summaries in one run.

library(plagueMetapop)
library(dplyr)
library(ggplot2)
library(here)
library(tidyr)

source(here::here("fadeout", "occupancy_functions.R"))
theme_set(theme_bw())

baseline <- list(
  R0 = 2.5, K = 1e4, r = 0.125, alpha = 1e-4, gamma = 1,
  n_patch = 200, dt = 0.1, t_max = 2000, I_outbreak = 10,
  initialization = "virgin_soil",
  initialization_seed = 101, simulation_seed = 102
)
parameter_values <- list(
  R0 = c(1.5, 2, 2.5, 3, 4, 5),
  K = c(1000, 3000, 10000, 30000, 100000),
  alpha = c(1e-5, 3e-5, 1e-4, 3e-4, 1e-3)
)
tau <- 50

make_scenarios <- function(baseline, parameter_values) {
  unlist(lapply(names(parameter_values), function(parameter) {
    lapply(parameter_values[[parameter]], function(value) {
      params <- baseline
      params[[parameter]] <- value
      list(varied_parameter = parameter, parameter_value = value, params = params)
    })
  }), recursive = FALSE)
}

format_parameter <- function(name, value) {
  if (name == "K") {
    sprintf("%s = %s", name, format(value, big.mark = ",", scientific = FALSE))
  } else {
    sprintf("%s = %g", name, value)
  }
}

fixed_parameter_subtitle <- function(varied_parameter, params) {
  fixed <- setdiff(c("R0", "K", "r", "alpha"), varied_parameter)
  paste0(
    paste(vapply(
      fixed, function(name) format_parameter(name, params[[name]]), character(1)
    ), collapse = "; "),
    sprintf(
      "; gamma = %g; n_patch = %d; dt = %g; t_max = %g; tau = %g\n",
      params$gamma, params$n_patch, params$dt, params$t_max, tau
    ),
    sprintf(
      "initialization = %s; S(0) = K - %g; I(0) = %g per patch",
      params$initialization, params$I_outbreak, params$I_outbreak
    )
  )
}

safe_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

integrate_trapezoid <- function(time, value) {
  if (length(time) < 2L) return(0)
  sum(diff(time) * (head(value, -1L) + tail(value, -1L)) / 2)
}

outdir <- here::here("fadeout", "output", "stochastic_episode_occupancy")
datadir <- file.path(outdir, "data")
figdir <- file.path(outdir, "figures")
dir.create(datadir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

plot_only <- "--plot-only" %in% commandArgs(trailingOnly = TRUE)

if (plot_only) {
  curves <- readRDS(file.path(
    datadir, "transient_source_pressure_timeseries.rds"
  ))
  summary_data <- read.csv(file.path(
    datadir, "transient_source_pressure_summary.csv"
  ))
  message("Plot-only mode: using existing episode-occupancy data")
} else {
  scenarios <- make_scenarios(baseline, parameter_values)
  gen <- compile_odin("euler_odin_def.R")
  results <- vector("list", length(scenarios))

  for (i in seq_along(scenarios)) {
    scenario <- scenarios[[i]]
    message(sprintf(
      "Running transient source-pressure audit: %s = %g (%d/%d)",
      scenario$varied_parameter, scenario$parameter_value, i, length(scenarios)
    ))
    simulation <- run_fadeout_occupancy(
      scenario$params, gen = gen, keep_patch_data = TRUE
    )
    contribution <- summarize_transient_source_pressure(
      simulation$infected,
      n_patch = scenario$params$n_patch,
      tau = tau,
      dt = scenario$params$dt
    )
    established <- summarize_established_occupancy(
      simulation$infected,
      n_patch = scenario$params$n_patch,
      tau = tau,
      dt = scenario$params$dt
    )
    results[[i]] <- list(
      varied_parameter = scenario$varied_parameter,
      parameter_value = scenario$parameter_value,
      params = scenario$params,
      tau = tau,
      contribution = contribution,
      occupancy_summary = established
    )
    rm(simulation, contribution, established)
    invisible(gc())
  }

  curves <- bind_rows(lapply(seq_along(results), function(i) {
    x <- results[[i]]
    x$contribution |>
      rename(time = step) |>
      mutate(
        result_id = i,
        varied_parameter = x$varied_parameter,
        parameter_value = x$parameter_value,
        R0 = x$params$R0, K = x$params$K, r = x$params$r,
        alpha = x$params$alpha, gamma = x$params$gamma, tau = x$tau,
        .before = 1
      )
  }))

  summary_data <- bind_rows(lapply(seq_along(results), function(i) {
    x <- results[[i]]
    z <- x$contribution
    ## Post-burnout time zero is the minimum retrospective persistent occupancy.
    minimum_index <- which.min(z$persistent_occupancy_fraction)
    minimum_time <- z$step[minimum_index]
    first500 <- z$step >= minimum_time & z$step <= minimum_time + 500
    integrated_transient_I <- integrate_trapezoid(z$step, z$transient_I)
    integrated_persistent_I <- integrate_trapezoid(z$step, z$persistent_I)
    cumulative_transient_pressure <-
      x$params$alpha * integrated_transient_I / x$params$n_patch
    cumulative_persistent_pressure <-
      x$params$alpha * integrated_persistent_I / x$params$n_patch
    cumulative_pressure_total <-
      cumulative_transient_pressure + cumulative_persistent_pressure
    integrated_classifiable_I <-
      integrated_transient_I + integrated_persistent_I
    cumulative_share <- if (cumulative_pressure_total > 0) {
      cumulative_transient_pressure / cumulative_pressure_total
    } else {
      NA_real_
    }
    integrated_I_share <- if (integrated_classifiable_I > 0) {
      integrated_transient_I / integrated_classifiable_I
    } else {
      NA_real_
    }
    data.frame(
      result_id = i,
      varied_parameter = x$varied_parameter,
      parameter_value = x$parameter_value,
      R0 = x$params$R0, K = x$params$K, r = x$params$r,
      alpha = x$params$alpha,
      tau = x$tau,
      occupancy_minimum_time = minimum_time,
      peak_transient_source_share = max(z$transient_source_share, na.rm = TRUE),
      transient_source_share_at_minimum =
        z$transient_source_share[minimum_index],
      mean_transient_source_share_first500 = safe_mean(
        z$transient_source_share[first500]
      ),
      mean_transient_source_share_total = safe_mean(z$transient_source_share),
      maximum_transient_I = max(z$transient_I),
      maximum_persistent_I = max(z$persistent_I),
      maximum_censored_I = max(z$censored_I),
      ratio_max_transient_to_persistent = if (max(z$persistent_I) > 0) {
        max(z$transient_I) / max(z$persistent_I)
      } else {
        NA_real_
      },
      integrated_transient_I = integrated_transient_I,
      integrated_persistent_I = integrated_persistent_I,
      cumulative_transient_pressure = cumulative_transient_pressure,
      cumulative_persistent_pressure = cumulative_persistent_pressure,
      cumulative_transient_source_share = cumulative_share,
      integrated_I_transient_share = integrated_I_share,
      cumulative_share_identity_error = cumulative_share - integrated_I_share
    )
  }))

  write.csv(
    curves,
    file.path(datadir, "transient_source_pressure_timeseries.csv"),
    row.names = FALSE
  )
  saveRDS(curves, file.path(
    datadir, "transient_source_pressure_timeseries.rds"
  ))
  write.csv(
    summary_data,
    file.path(datadir, "transient_source_pressure_summary.csv"),
    row.names = FALSE
  )

  ## Save the established target in the same list structure expected by the
  ## analytical comparison scripts.
  established_results <- lapply(results, function(x) {
    list(
      varied_parameter = x$varied_parameter,
      parameter_value = x$parameter_value,
      params = x$params,
      occupancy_summary = x$occupancy_summary
    )
  })
  established_curves <- bind_rows(lapply(established_results, function(x) {
    x$occupancy_summary |>
      mutate(
        varied_parameter = x$varied_parameter,
        parameter_value = x$parameter_value,
        R0 = x$params$R0, K = x$params$K, r = x$params$r,
        alpha = x$params$alpha, .before = 1
      )
  }))
  saveRDS(established_results,
    file.path(datadir, "established_occupancy_results.rds"))
  write.csv(established_curves,
    file.path(datadir, "established_occupancy_curves.csv"), row.names = FALSE)
}

special_R0 <- curves |>
  filter(
    varied_parameter == "R0",
    parameter_value %in% c(2.5, 3),
    time <= max(time) - tau
  ) |>
  mutate(R0_label = factor(
    vapply(parameter_value, function(value) format_parameter("R0", value),
           character(1)),
    levels = c("R0 = 2.5", "R0 = 3")
  ))
p_special_R0 <- ggplot(
  special_R0,
  aes(time, transient_source_share, colour = R0_label)
) +
  geom_line(linewidth = 0.7, na.rm = TRUE) +
  scale_colour_manual(values = c("#0072B2", "#D55E00")) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(
    x = "Time (disease generations)",
    y = "Transient source share of colonization pressure",
    colour = NULL,
    title = "Transient source pressure: baseline R0 versus R0 = 3",
    subtitle = fixed_parameter_subtitle("R0", baseline)
  ) +
  theme(legend.position = "right", plot.subtitle = element_text(size = 9))
ggsave(
  file.path(figdir, "transient_source_share_R0_2p5_vs_3.pdf"),
  p_special_R0, width = 10, height = 6
)

special_components <- special_R0 |>
  select(time, R0_label, persistent_I, transient_I) |>
  pivot_longer(
    c(persistent_I, transient_I),
    names_to = "component", values_to = "infected_hosts"
  ) |>
  mutate(component = factor(
    component,
    levels = c("persistent_I", "transient_I"),
    labels = c("Persistent-episode I", "Transient-episode I")
  ))
p_special_components <- ggplot(
  special_components,
  aes(time, infected_hosts, fill = component)
) +
  geom_area(position = "stack", alpha = 0.75) +
  geom_line(
    data = special_R0,
    aes(time, total_I, group = R0_label),
    inherit.aes = FALSE, colour = "black", linewidth = 0.5
  ) +
  facet_wrap(~R0_label, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("#0072B2", "#D55E00")) +
  scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
  labs(
    x = "Time (disease generations)", y = "Infected hosts", fill = NULL,
    title = "Persistent and transient source components: R0 = 2.5 versus 3",
    subtitle = fixed_parameter_subtitle("R0", baseline),
    caption = "Stacked areas equal total_I; only t = 0--1950 is shown."
  ) +
  theme(
    legend.position = "top", plot.subtitle = element_text(size = 9),
    plot.caption = element_text(hjust = 0, size = 8)
  )
ggsave(
  file.path(figdir, "infected_components_R0_2p5_vs_3.pdf"),
  p_special_components, width = 10, height = 8
)

message("Finished. New outputs saved in: ", outdir)
