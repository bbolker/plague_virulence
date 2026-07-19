## Primary comparison of forward-window established stochastic occupancy with
## the unchanged analytical logistic approximation. This script reads existing
## established-occupancy summaries and does not run metapopulation simulations.

library(dplyr)
library(ggplot2)
library(here)
library(mgcv)
library(plagueMetapop)

theme_set(theme_bw())

input_file <- here::here(
  "fadeout", "output", "patch_occupancy_established", "data",
  "established_occupancy_results.rds"
)
single_patch_file <- here::here(
  "odin", "sharcnet", "outputs",
  "euler_onepatch_onestrain_extinct_logistic_continuous.rds"
)
outdir <- here::here(
  "fadeout", "output", "patch_occupancy_established",
  "analytical_comparison"
)
datadir <- file.path(outdir, "data")
figdir <- file.path(outdir, "figures")

for (filename in c(input_file, single_patch_file)) {
  if (!file.exists(filename)) stop("Required input not found: ", filename)
}
dir.create(datadir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

results <- readRDS(input_file)
single_patch <- readRDS(single_patch_file)
required_single <- c("R0", "K", "ext_prob.I1")
if (!is.data.frame(single_patch) || nrow(single_patch) == 0L ||
    length(setdiff(required_single, names(single_patch))) > 0L) {
  stop("Single-patch data must be non-empty and contain R0, K, ext_prob.I1")
}
if (!is.list(results) || length(results) == 0L) {
  stop("Established-occupancy input must be a non-empty list")
}

## Match compare_patch_occupancy_approximation.R exactly: the GAM predicts
## extinction probability, P1 is one minus that prediction, and numerical
## evaluation is bounded away from zero and one.
gam_fit <- mgcv::gam(
  ext_prob.I1 ~ te(R0, K, k = c(12, 12)),
  data = single_patch,
  method = "REML"
)
probability_tolerance <- 1e-6
## Use one conservative persistence definition. With dt=0.1 in the current
## scenarios, tau=50 corresponds to 500 stored simulation steps.
tau_values <- 50

format_parameter <- function(name, value) {
  if (name == "K") {
    sprintf("%s = %s", name, format(value, big.mark = ",", scientific = FALSE))
  } else {
    sprintf("%s = %g", name, value)
  }
}

fixed_parameter_subtitle <- function(varied_parameter, params) {
  fixed <- setdiff(c("R0", "K", "r", "alpha"), varied_parameter)
  biological <- paste(vapply(
    fixed, function(name) format_parameter(name, params[[name]]), character(1)
  ), collapse = "; ")
  initialization <- if (identical(params$initialization, "virgin_soil")) {
    sprintf(
      "initialization = virgin_soil; S(0) = K - %g; I(0) = %g per patch",
      params$I_outbreak, params$I_outbreak
    )
  } else {
    sprintf("initialization = %s", params$initialization)
  }
  paste0(
    biological,
    sprintf(
      "; gamma = %g; n_patch = %d; dt = %g; t_max = %g\n",
      params$gamma, params$n_patch, params$dt, params$t_max
    ),
    initialization
  )
}

first_time_at_or_above <- function(time, value, level) {
  index <- which(!is.na(value) & value >= level)[1L]
  if (is.na(index)) NA_real_ else time[index]
}

analyse_one_tau <- function(result, result_id, tau_value) {
  required_result <- c(
    "varied_parameter", "parameter_value", "params", "occupancy_summary"
  )
  missing_result <- setdiff(required_result, names(result))
  if (length(missing_result) > 0L) {
    stop("Result ", result_id, " is missing fields: ",
         paste(missing_result, collapse = ", "))
  }

  required_summary <- c("step", "tau", "occupancy_fraction_established")
  missing_summary <- setdiff(required_summary, names(result$occupancy_summary))
  if (length(missing_summary) > 0L) {
    stop("Result ", result_id, " occupancy summary is missing columns: ",
         paste(missing_summary, collapse = ", "))
  }

  simulation <- result$occupancy_summary |>
    filter(.data$tau == .env$tau_value) |>
    select(step, occupancy_established = occupancy_fraction_established) |>
    filter(!is.na(occupancy_established))
  if (nrow(simulation) == 0L) {
    stop("No nonmissing established occupancy for result ", result_id,
         ", tau=", tau_value)
  }

  params <- result$params
  extinction_prediction <- as.numeric(predict(
    gam_fit,
    newdata = data.frame(R0 = params$R0, K = params$K),
    type = "response"
  ))
  P1_GAM_raw <- 1 - extinction_prediction
  P1 <- pmin(pmax(P1_GAM_raw, probability_tolerance),
             1 - probability_tolerance)
  exact_row <- single_patch |>
    filter(R0 == params$R0, K == params$K)
  P1_exact_grid <- if (nrow(exact_row) == 1L) {
    1 - exact_row$ext_prob.I1
  } else {
    NA_real_
  }

  equilibrium <- plagueMetapop::ode_eq(
    beta = params$R0,
    gamma = params$gamma,
    K = params$K,
    r = params$r,
    logistic_growth = 1
  )
  I_star <- unname(equilibrium["eq_I"])
  if (!is.finite(I_star) || I_star <= 0) {
    stop("Invalid endemic I_star for result ", result_id, ": ", I_star)
  }
  lambda <- params$alpha * I_star * P1

  ## Alignment uses the minimum of this established curve, never raw
  ## occupancy. The analytical initial condition remains p(0)=P1 and is not
  ## fitted to the observed minimum.
  minimum_index <- which.min(simulation$occupancy_established)
  time_shift <- simulation$step[minimum_index]
  comparison <- simulation |>
    filter(step >= time_shift) |>
    transmute(
      time_post_burnout = step - time_shift,
      occupancy_established,
      occupancy_approximation = 1 / (
        1 + ((1 - P1) / P1) * exp(-lambda * time_post_burnout)
      )
    ) |>
    mutate(residual = occupancy_established - occupancy_approximation)
  early <- comparison |>
    filter(time_post_burnout <= 500)

  parameter_label <- format_parameter(
    result$varied_parameter, result$parameter_value
  )
  curves <- bind_rows(
    comparison |>
      transmute(
        time_post_burnout,
        occupancy = occupancy_established,
        trajectory = "Established stochastic occupancy"
      ),
    comparison |>
      transmute(
        time_post_burnout,
        occupancy = occupancy_approximation,
        trajectory = "Analytical approximation"
      )
  ) |>
    mutate(
      result_id = result_id,
      varied_parameter = result$varied_parameter,
      parameter_value = result$parameter_value,
      parameter_label = parameter_label,
      tau = tau_value
    )

  correlation <- if (sd(comparison$occupancy_established) == 0 ||
                       sd(comparison$occupancy_approximation) == 0) {
    NA_real_
  } else {
    cor(comparison$occupancy_established,
        comparison$occupancy_approximation)
  }
  diagnostics <- data.frame(
    varied_parameter = result$varied_parameter,
    parameter_value = result$parameter_value,
    R0 = params$R0,
    K = params$K,
    r = params$r,
    gamma = params$gamma,
    alpha = params$alpha,
    tau = tau_value,
    P1_GAM_raw = P1_GAM_raw,
    P1_used = P1,
    P1_exact_grid = P1_exact_grid,
    I_star = I_star,
    lambda = lambda,
    simulation_minimum_established =
      simulation$occupancy_established[minimum_index],
    time_shift_established = time_shift,
    RMSE = sqrt(mean(comparison$residual^2)),
    RMSE_first_500 = sqrt(mean(early$residual^2)),
    max_absolute_deviation = max(abs(comparison$residual)),
    max_absolute_deviation_first_500 = max(abs(early$residual)),
    correlation = correlation,
    simulation_time_to_90pct = first_time_at_or_above(
      comparison$time_post_burnout, comparison$occupancy_established, 0.90
    ),
    approximation_time_to_90pct = first_time_at_or_above(
      comparison$time_post_burnout, comparison$occupancy_approximation, 0.90
    ),
    simulation_time_to_95pct = first_time_at_or_above(
      comparison$time_post_burnout, comparison$occupancy_established, 0.95
    ),
    approximation_time_to_95pct = first_time_at_or_above(
      comparison$time_post_burnout, comparison$occupancy_approximation, 0.95
    ),
    number_of_nonmissing_time_points = nrow(simulation),
    stringsAsFactors = FALSE
  )

  list(curves = curves, diagnostics = diagnostics, params = params)
}

analyses <- unlist(lapply(seq_along(results), function(i) {
  lapply(tau_values, function(tau_value) {
    analyse_one_tau(results[[i]], i, tau_value)
  })
}), recursive = FALSE)
all_curves <- bind_rows(lapply(analyses, `[[`, "curves"))
all_diagnostics <- bind_rows(lapply(analyses, `[[`, "diagnostics"))

palette <- c(
  "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"
)

for (varied_parameter in unique(all_curves$varied_parameter)) {
  parameter_name <- varied_parameter
  group_rows <- which(all_diagnostics$varied_parameter == parameter_name)
  group_diagnostics <- all_diagnostics[group_rows, , drop = FALSE]
  label_levels <- unique(group_diagnostics[
    order(group_diagnostics$parameter_value), c("parameter_value")
  ])
  label_levels <- vapply(label_levels, function(value) {
    format_parameter(parameter_name, value)
  }, character(1))
  fixed_params <- analyses[[group_rows[1]]]$params

  for (tau_value in tau_values) {
    plot_data <- all_curves |>
      filter(.data$varied_parameter == .env$parameter_name,
             .data$tau == .env$tau_value) |>
      mutate(parameter_label = factor(parameter_label, levels = label_levels))
    p <- ggplot(
      plot_data,
      aes(
        time_post_burnout, occupancy,
        colour = parameter_label,
        linetype = trajectory,
        group = interaction(parameter_label, trajectory)
      )
    ) +
      geom_line(linewidth = 0.65) +
      scale_colour_manual(values = palette[seq_along(label_levels)]) +
      scale_linetype_manual(values = c(
        "Established stochastic occupancy" = "solid",
        "Analytical approximation" = "22"
      )) +
      scale_y_continuous(
        limits = c(0, 1), breaks = seq(0, 1, 0.2),
        labels = scales::percent_format(accuracy = 1)
      ) +
      labs(
        x = "Time since established-occupancy minimum",
        y = "Fraction of established patches",
        colour = sprintf("Varied parameter: %s", parameter_name),
        linetype = NULL,
        title = sprintf(
          "Established stochastic occupancy vs analytical approximation: %s",
          parameter_name
        ),
        subtitle = paste0(
          "tau = ", tau_value, "; ",
          fixed_parameter_subtitle(parameter_name, fixed_params)
        ),
        caption = paste0(
          "Time zero is the minimum nonmissing established occupancy for this tau. ",
          "The analytical curve starts at p(0) = P1; P1 is not refitted."
        )
      ) +
      theme(
        legend.position = "right",
        plot.subtitle = element_text(size = 9),
        plot.caption = element_text(hjust = 0, size = 8)
      )
    ggsave(
      file.path(figdir, sprintf(
        "analytical_vs_established_%s_tau_%g.pdf",
        parameter_name, tau_value
      )),
      p, width = 11, height = 6.5
    )
  }
}

write.csv(
  all_diagnostics,
  file.path(datadir, "established_occupancy_approximation_diagnostics.csv"),
  row.names = FALSE
)
write.csv(
  all_curves,
  file.path(datadir, "established_occupancy_approximation_curves.csv"),
  row.names = FALSE
)

cat("Primary comparison: established stochastic occupancy vs analytical approximation\n")
cat("Scenario-tau comparisons: ", nrow(all_diagnostics), "\n", sep = "")
print(all_diagnostics |>
  group_by(tau) |>
  summarise(
    median_RMSE = median(RMSE),
    median_RMSE_first_500 = median(RMSE_first_500),
    median_abs_P1_minus_minimum = median(abs(P1_used -
      simulation_minimum_established)),
    .groups = "drop"
  ))
cat("Outputs: ", outdir, "\n", sep = "")
