## Compare the provisional half-period transient duration with a hybrid
## deterministic endpoint: use the post-peak downward I=1 crossing when the
## first trough is <= 1, otherwise use the first trough itself.
##
## Run from the repository root:
##   Rscript fadeout/compare_transient_duration_methods.R

library(dplyr)
library(ggplot2)
library(here)
library(cowplot)
library(tidyr)

source(here::here("fadeout", "two_state_functions.R"))
theme_set(theme_bw())

baseline <- list(R0 = 2.5, K = 10000, r = 0.125, I0 = 10)
parameter_values <- list(
  R0 = c(1.5, 2, 2.5, 3, 4, 5),
  K = c(1000, 3000, 10000, 30000, 100000),
  r = c(0.05, 0.1, 0.125, 0.2, 0.4)
)

outdir <- here::here("fadeout", "output", "transient_duration_methods")
datadir <- file.path(outdir, "data")
figdir <- file.path(outdir, "figures")
dir.create(datadir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

run_pair <- function(params, integration_dt = 0.01) {
  old <- compute_transient_outbreak_summary(
    R0 = params$R0, K = params$K, r = params$r, gamma = 1,
    I0 = params$I0, c = 0.5, integration_dt = integration_dt
  )
  new <- compute_transient_outbreak_summary_I1(
    R0 = params$R0, K = params$K, r = params$r, I0 = params$I0,
    integration_dt = integration_dt, initial_tmax = 50,
    maximum_tmax = 800
  )
  list(old = old, new = new)
}

baseline_result <- run_pair(baseline)
old <- baseline_result$old
new <- baseline_result$new

baseline_table <- bind_rows(
  data.frame(
    method = "old_half_Tosc",
    endpoint_type = "half_Tosc",
    T_T = old$T_T,
    integral_I = old$integral_I,
    Ibar_T = old$Ibar_T,
    I_peak = old$I_peak,
    time_of_I_peak = old$time_of_I_peak,
    endpoint_found = TRUE,
    time_of_first_trough = new$time_of_first_trough,
    I_first_trough = new$I_first_trough
  ),
  data.frame(
    method = "I1_or_first_trough",
    endpoint_type = new$transient_endpoint,
    T_T = new$T_T,
    integral_I = new$integral_I,
    Ibar_T = new$Ibar_T,
    I_peak = new$I_peak,
    time_of_I_peak = new$time_of_I_peak,
    endpoint_found = new$endpoint_found,
    time_of_first_trough = new$time_of_first_trough,
    I_first_trough = new$I_first_trough
  )
)
write.csv(
  baseline_table,
  file.path(datadir, "transient_duration_baseline_comparison.csv"),
  row.names = FALSE
)

## Baseline trajectory: retain only the first cycle and a small margin after
## its trough. The full adaptively extended trajectory remains available from
## the function return object but is not needed for this diagnostic figure.
plot_end <- if (new$trough_found) {
  min(max(new$time_of_first_trough * 1.25, old$T_T * 1.25), 100)
} else {
  min(new$integration_horizon, 100)
}
trajectory_plot <- new$trajectory |>
  filter(time <= plot_end, is.finite(I), I > 0)

p_baseline <- ggplot(trajectory_plot, aes(time, I)) +
  geom_line(linewidth = 0.75, colour = "#0072B2") +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey35") +
  geom_vline(
    xintercept = old$T_T, colour = "#D55E00",
    linewidth = 0.75, linetype = "dotdash"
  ) +
  geom_point(
    data = data.frame(time = new$time_of_I_peak, I = new$I_peak),
    colour = "#009E73", size = 2.8
  ) +
  geom_point(
    data = data.frame(
      time = new$time_of_first_trough, I = new$I_first_trough
    ),
    colour = "#CC79A7", size = 2.8
  ) +
  geom_vline(
    xintercept = new$T_T, colour = "black", linewidth = 0.9
  ) +
  geom_point(
    data = data.frame(time = new$T_T, I = new$I_at_T_T),
    colour = "black", size = 3
  ) +
  scale_y_log10(labels = scales::label_number(big.mark = ",")) +
  labs(
    x = "Time (disease generations)",
    y = "Deterministic infected hosts I(t) (log scale)",
    title = "Transient-duration definitions at the baseline parameters",
    subtitle = paste0(
      "R0 = 2.5; K = 10,000; r = 0.125; I(0) = 10; ",
      "logistic susceptible growth"
    ),
    caption = paste0(
      "Green: first peak; magenta: first trough; orange dot-dash: old ",
      "T_T = 0.5 T_osc; dashed horizontal: I = 1.\n",
      "Black: selected hybrid endpoint (", new$transient_endpoint,
      "); T_T = ", format(new$T_T, digits = 5), "."
    )
  ) +
  theme(
    plot.subtitle = element_text(size = 9),
    plot.caption = element_text(hjust = 0, size = 8)
  )
ggsave(
  file.path(figdir, "transient_duration_baseline_trajectory.pdf"),
  p_baseline, width = 10, height = 7
)

make_scenarios <- function() {
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

scenarios <- make_scenarios()
ofat_results <- bind_rows(lapply(seq_along(scenarios), function(i) {
  scenario <- scenarios[[i]]
  message(sprintf(
    "Computing %s = %g (%d/%d)",
    scenario$varied_parameter, scenario$parameter_value,
    i, length(scenarios)
  ))
  result <- run_pair(scenario$params)
  data.frame(
    varied_parameter = scenario$varied_parameter,
    parameter_value = scenario$parameter_value,
    R0 = scenario$params$R0,
    K = scenario$params$K,
    r = scenario$params$r,
    I0 = scenario$params$I0,
    old_T_osc = result$old$T_osc,
    old_T_T = result$old$T_T,
    old_integral_I = result$old$integral_I,
    old_Ibar_T = result$old$Ibar_T,
    endpoint_found = result$new$endpoint_found,
    endpoint_type = result$new$transient_endpoint,
    hybrid_status = result$new$status,
    hybrid_T_T = result$new$T_T,
    hybrid_integral_I = result$new$integral_I,
    hybrid_Ibar_T = result$new$Ibar_T,
    I_peak = result$new$I_peak,
    time_of_I_peak = result$new$time_of_I_peak,
    time_of_first_trough = result$new$time_of_first_trough,
    I_first_trough = result$new$I_first_trough,
    I_at_hybrid_T_T = result$new$I_at_T_T,
    integration_horizon = result$new$integration_horizon
  )
}))
write.csv(
  ofat_results,
  file.path(datadir, "transient_duration_ofat_comparison.csv"),
  row.names = FALSE
)

format_parameter <- function(name, value) {
  if (name == "K") {
    sprintf("K = %s", format(value, scientific = FALSE, big.mark = ","))
  } else {
    sprintf("%s = %g", name, value)
  }
}

fixed_parameter_text <- function(varied_parameter) {
  fixed <- setdiff(c("R0", "K", "r"), varied_parameter)
  paste(vapply(
    fixed,
    function(name) format_parameter(name, baseline[[name]]),
    character(1)
  ), collapse = "; ")
}

make_method_figure <- function(parameter) {
  data <- ofat_results |>
    filter(varied_parameter == parameter) |>
    arrange(parameter_value)

  duration <- data |>
    select(parameter_value, old_T_T, hybrid_T_T) |>
    pivot_longer(
      c(old_T_T, hybrid_T_T), names_to = "method", values_to = "value"
    ) |>
    mutate(method = recode(
      method,
      old_T_T = "Old: 0.5 T_osc",
      hybrid_T_T = "Hybrid: I=1 or first trough"
    ))

  infected_load <- data |>
    select(parameter_value, old_Ibar_T, hybrid_Ibar_T) |>
    pivot_longer(
      c(old_Ibar_T, hybrid_Ibar_T),
      names_to = "method", values_to = "value"
    ) |>
    mutate(method = recode(
      method,
      old_Ibar_T = "Old: 0.5 T_osc",
      hybrid_Ibar_T = "Hybrid: I=1 or first trough"
    ))

  colours <- c(
    "Old: 0.5 T_osc" = "#D55E00",
    "Hybrid: I=1 or first trough" = "#0072B2"
  )

  p_duration <- ggplot(
    duration |> filter(is.finite(value)),
    aes(parameter_value, value, colour = method)
  ) +
    geom_line(aes(group = method), linewidth = 0.7) +
    geom_point(size = 2.1) +
    scale_colour_manual(values = colours) +
    labs(
      x = NULL, y = expression("Transient duration " * T[T]),
      colour = NULL,
      title = sprintf(
        "Transient-duration methods: varying %s", parameter
      ),
      subtitle = paste0(
        fixed_parameter_text(parameter), "; I(0) = 10; ",
        "time in disease generations"
      )
    ) +
    theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      legend.position = "top",
      plot.subtitle = element_text(size = 9),
      plot.margin = margin(8, 8, 4, 80, unit = "pt")
    )

  p_Ibar <- ggplot(
    infected_load |> filter(is.finite(value)),
    aes(parameter_value, value, colour = method)
  ) +
    geom_line(aes(group = method), linewidth = 0.7) +
    geom_point(size = 2.1) +
    scale_colour_manual(values = colours) +
    labs(
      x = NULL, y = expression("Mean transient infected load " * bar(I)[T]),
      colour = NULL
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      plot.margin = margin(4, 8, 4, 80, unit = "pt")
    )

  p_trough <- ggplot(data, aes(parameter_value, I_first_trough)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_line(linewidth = 0.7, colour = "#CC79A7") +
    geom_point(size = 2.1, colour = "#CC79A7") +
    scale_y_log10() +
    labs(
      x = NULL, y = "First-trough I (log scale)"
    ) +
    theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      plot.margin = margin(4, 8, 4, 80, unit = "pt")
    )

  endpoint_levels <- c("I1_crossing", "first_trough", "not_found")
  endpoint_data <- data |>
    mutate(
      endpoint_plot = ifelse(endpoint_found, endpoint_type, "not_found"),
      endpoint_plot = factor(endpoint_plot, levels = endpoint_levels)
    )
  p_endpoint <- ggplot(
    endpoint_data, aes(parameter_value, endpoint_plot, colour = endpoint_plot)
  ) +
    geom_point(size = 3) +
    scale_colour_manual(
      values = c(
        I1_crossing = "#009E73", first_trough = "#CC79A7",
        not_found = "red3"
      ),
      drop = FALSE
    ) +
    labs(
      x = parameter, y = "Selected endpoint", colour = NULL,
      caption = "Rule: first trough I <= 1 uses the downward I=1 crossing; otherwise uses the first trough."
    ) +
    theme(
      legend.position = "none",
      plot.caption = element_text(hjust = 0, size = 8),
      plot.margin = margin(4, 8, 8, 80, unit = "pt")
    )

  if (parameter == "K") {
    K_scale <- scale_x_log10(
      breaks = c(1000, 3000, 10000, 30000, 100000),
      labels = scales::label_number(big.mark = ",")
    )
    p_duration <- p_duration + K_scale
    p_Ibar <- p_Ibar + K_scale
    p_trough <- p_trough + K_scale
    p_endpoint <- p_endpoint + K_scale + labs(x = "K (log scale)")
  }

  plot_grid(
    p_duration, p_Ibar, p_trough, p_endpoint,
    ncol = 1, rel_heights = c(1, 0.85, 0.85, 0.95),
    align = "v", axis = "l"
  )
}

for (parameter in names(parameter_values)) {
  ggsave(
    file.path(
      figdir, sprintf("transient_duration_compare_%s.pdf", parameter)
    ),
    make_method_figure(parameter), width = 10, height = 12
  )
}

## Lightweight validations.
old_expected <- c(
  T_osc = 15.558639889, T_T = 7.7793199445,
  integral_I = 8863.5128, Ibar_T = 1139.3690
)
old_unchanged <- isTRUE(all.equal(
  unlist(old[names(old_expected)]), old_expected,
  tolerance = 2e-6
))

fine_baseline <- compute_transient_outbreak_summary_I1(
  2.5, 10000, 0.125, I0 = 10, integration_dt = 0.005,
  initial_tmax = 50, maximum_tmax = 800
)
crossing_example <- compute_transient_outbreak_summary_I1(
  2.5, 1000, 0.125, I0 = 10, integration_dt = 0.01,
  initial_tmax = 50, maximum_tmax = 800
)
crossing_example_fine <- compute_transient_outbreak_summary_I1(
  2.5, 1000, 0.125, I0 = 10, integration_dt = 0.005,
  initial_tmax = 50, maximum_tmax = 800
)

crossing_checks <- c(
  baseline_endpoint_found = new$endpoint_found,
  baseline_uses_trough = identical(new$transient_endpoint, "first_trough"),
  baseline_trough_above_one = new$I_first_trough > 1,
  crossing_endpoint_found = crossing_example$endpoint_found,
  crossing_uses_I1 = identical(
    crossing_example$transient_endpoint, "I1_crossing"
  ),
  crossing_trough_at_or_below_one = crossing_example$I_first_trough <= 1,
  crossing_after_peak =
    crossing_example$T_T > crossing_example$time_of_I_peak,
  before_above_one = crossing_example$I_before_crossing > 1,
  after_at_or_below_one = crossing_example$I_after_crossing <= 1,
  interpolated_at_one = abs(crossing_example$I_at_T_T - 1) < 1e-8,
  positive_Ibar = crossing_example$Ibar_T > 0,
  Ibar_below_peak =
    crossing_example$Ibar_T <= crossing_example$I_peak + 1e-7,
  stable_T_T =
    abs(crossing_example$T_T - crossing_example_fine$T_T) < 0.01,
  stable_Ibar =
    abs(crossing_example$Ibar_T - crossing_example_fine$Ibar_T) /
      crossing_example_fine$Ibar_T < 1e-3,
  stable_baseline_trough =
    abs(new$T_T - fine_baseline$T_T) < 0.01 &&
      abs(new$Ibar_T - fine_baseline$Ibar_T) /
        fine_baseline$Ibar_T < 1e-3,
  peak_S_crossing =
    new$peak_S_left > new$equilibrium_S &&
      new$peak_S_right <= new$equilibrium_S,
  trough_S_crossing =
    new$trough_S_left < new$equilibrium_S &&
      new$trough_S_right >= new$equilibrium_S,
  interpolated_peak_S =
    abs(new$S_at_peak - new$equilibrium_S) < 1e-7,
  interpolated_trough_S =
    abs(new$S_at_first_trough - new$equilibrium_S) < 1e-7,
  normalized_time_units =
    identical(crossing_example$time_units, "disease_generations"),
  old_function_unchanged = old_unchanged
)
if (!all(crossing_checks)) {
  stop("Transient-duration validation failed for: ",
       paste(names(crossing_checks)[!crossing_checks], collapse = ", "))
}

endpoint_summary <- ofat_results |>
  count(endpoint_type, endpoint_found, name = "n_scenarios")
write.csv(
  endpoint_summary,
  file.path(datadir, "transient_duration_endpoint_summary.csv"),
  row.names = FALSE
)

cat("Validation checks passed: ", length(crossing_checks), "/",
    length(crossing_checks), "\n", sep = "")
cat("Baseline old T_T: ", old$T_T, "\n", sep = "")
cat("Baseline hybrid endpoint: ", new$transient_endpoint,
    "; hybrid T_T: ", ifelse(is.na(new$T_T), "NA", new$T_T),
    "; first trough I: ", new$I_first_trough, "\n", sep = "")
cat("OFAT endpoint counts:\n")
print(endpoint_summary)
cat("Outputs saved under: ", outdir, "\n", sep = "")
