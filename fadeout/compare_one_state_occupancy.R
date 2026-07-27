## Compare the analytical patch-occupancy approximation with every existing
## trajectory in the one-factor-at-a-time stochastic metapopulation analysis.
## This script reads existing outputs only; it does not run simulations.

library(dplyr)
library(ggplot2)
library(here)
library(mgcv)
library(plagueMetapop)

source(here::here("fadeout", "occupancy_functions.R"))
theme_set(theme_bw())

single_patch_file <- here::here(
  "odin", "sharcnet", "outputs",
  "euler_onepatch_onestrain_extinct_logistic_continuous_demoggrid.rds"
)
metapop_file <- here::here(
  "fadeout", "output", "stochastic_patch_occupancy", "data",
  "stochastic_patch_occupancy_results.rds"
)
established_file <- here::here(
  "fadeout", "output", "stochastic_episode_occupancy", "data",
  "established_occupancy_results.rds"
)
outdir <- here::here("fadeout", "output", "one_state_occupancy")
diagnostic_file <- file.path(outdir, "one_state_occupancy_diagnostics.csv")

for (filename in c(single_patch_file, metapop_file, established_file)) {
  if (!file.exists(filename)) stop("Required input file not found: ", filename)
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

single_patch <- readRDS(single_patch_file)
required_single <- c("R0", "K", "r", "ext_prob.I1")
missing_single <- setdiff(required_single, names(single_patch))
if (!is.data.frame(single_patch) || nrow(single_patch) == 0L) {
  stop("Single-patch input must be a non-empty data frame")
}
if (length(missing_single) > 0L) {
  stop("Single-patch data are missing columns: ",
       paste(missing_single, collapse = ", "))
}

metapop_results <- readRDS(metapop_file)
established_results <- readRDS(established_file)
if (!is.list(metapop_results) || length(metapop_results) == 0L) {
  stop("Metapopulation input must be a non-empty list")
}
if (!is.list(established_results) ||
    length(established_results) != length(metapop_results)) {
  stop("Established and raw metapopulation inputs must be equal-length lists")
}

## Exact R0 x K x r grid values are preferred. Non-grid R0/K combinations use
## the two-dimensional GAM fitted within the matching observed r slice.
estimate_P1 <- make_P1_demoggrid_estimator(single_patch)

probability_tolerance <- 1e-6

format_parameter <- function(name, value) {
  if (name == "K") {
    sprintf("%s = %s", name, format(value, scientific = FALSE, big.mark = ","))
  } else {
    sprintf("%s = %g", name, value)
  }
}

fixed_parameter_subtitle <- function(group, params) {
  fixed <- setdiff(c("R0", "K", "r", "alpha"), group)
  paste0(
    paste(vapply(
      fixed,
      function(name) format_parameter(name, params[[name]]),
      character(1)
    ), collapse = "; "),
    sprintf(
      "; gamma = %g; n_patch = %d; dt = %g; t_max = %g\n",
      params$gamma, params$n_patch, params$dt, params$t_max
    ),
    sprintf(
      "initialization = %s; S(0) = K - %g; I(0) = %g per patch; tau = 50",
      params$initialization, params$I_outbreak, params$I_outbreak
    )
  )
}

first_time_at_or_above <- function(time, value, level) {
  index <- which(value >= level)[1L]
  if (is.na(index)) NA_real_ else time[index]
}

analyse_one <- function(result, established, result_id) {
  required_result <- c(
    "varied_parameter", "parameter_value", "params", "raw_meta_summary"
  )
  missing_result <- setdiff(required_result, names(result))
  if (length(missing_result) > 0L) {
    stop("Metapopulation result ", result_id, " is missing fields: ",
         paste(missing_result, collapse = ", "))
  }

  if (!identical(result$varied_parameter, established$varied_parameter) ||
      !isTRUE(all.equal(result$parameter_value,
        established$parameter_value))) {
    stop("Raw/established scenario mismatch for result ", result_id)
  }
  params <- result$params
  raw_simulation <- result$raw_meta_summary
  required_meta <- c("step", "occupied_patches", "occupancy_fraction")
  missing_meta <- setdiff(required_meta, names(raw_simulation))
  if (length(missing_meta) > 0L) {
    stop("Trajectory ", result_id, " is missing columns: ",
         paste(missing_meta, collapse = ", "))
  }
  simulation <- established$occupancy_summary |>
    filter(tau == 50, !is.na(occupancy_fraction_established)) |>
    select(step, occupancy_fraction = occupancy_fraction_established)
  if (nrow(simulation) == 0L) {
    stop("No tau=50 established occupancy for result ", result_id)
  }

  P1_estimate <- estimate_P1(params$R0, params$K, params$r)
  P1_raw <- P1_estimate$P1_raw
  P1 <- pmin(pmax(P1_raw, probability_tolerance),
             1 - probability_tolerance)
  P1_was_bounded <- !isTRUE(all.equal(P1_raw, P1))

  equilibrium <- plagueMetapop::ode_eq(
    beta = params$R0,
    gamma = params$gamma,
    K = params$K,
    r = params$r,
    logistic_growth = 1
  )
  I_star <- unname(equilibrium["eq_I"])
  if (!is.finite(I_star) || I_star <= 0) {
    stop("Invalid endemic I* for result ", result_id, ": ", I_star)
  }

  ## Use tau=50 established occupancy as the stochastic comparison target.
  ## Retain the raw-occupancy minimum as the common post-burnout time origin
  ## used by the current two-state comparison.
  minimum_index <- which.min(raw_simulation$occupancy_fraction)
  time_shift <- raw_simulation$step[minimum_index]
  post_burnout <- simulation |>
    filter(step >= time_shift) |>
    transmute(
      time_post_burnout = step - time_shift,
      occupancy_simulation = occupancy_fraction
    )

  lambda <- params$alpha * I_star * P1
  comparison <- post_burnout |>
    mutate(
      occupancy_approximation = 1 / (
        1 + ((1 - P1) / P1) * exp(-lambda * time_post_burnout)
      ),
      residual = occupancy_simulation - occupancy_approximation
    )
  early <- comparison |>
    filter(time_post_burnout <= 500)

  parameter_label <- format_parameter(
    result$varied_parameter, result$parameter_value
  )
  curves <- bind_rows(
    comparison |>
      transmute(
        time_post_burnout,
        occupancy = occupancy_simulation,
        trajectory = "Established stochastic occupancy (tau = 50)"
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
      group = result$varied_parameter,
      parameter_value = result$parameter_value,
      parameter_label = parameter_label
    )

  diagnostics <- data.frame(
    result_id = result_id,
    group = result$varied_parameter,
    parameter_value = result$parameter_value,
    R0 = params$R0,
    K = params$K,
    r = params$r,
    gamma = params$gamma,
    alpha = params$alpha,
    tau = 50,
    P1_raw = P1_raw,
    P1_used = P1,
    P1_was_bounded = P1_was_bounded,
    P1_exact_grid = P1_estimate$P1_exact_grid,
    P1_source = P1_estimate$source,
    I_star = I_star,
    lambda = lambda,
    simulation_established_at_shift = post_burnout$occupancy_simulation[1],
    raw_simulation_minimum = raw_simulation$occupancy_fraction[minimum_index],
    time_shift = time_shift,
    n_simulation_replicates = 1L,
    RMSE = sqrt(mean(comparison$residual^2)),
    max_absolute_deviation = max(abs(comparison$residual)),
    correlation = if (sd(comparison$occupancy_simulation) == 0 ||
                      sd(comparison$occupancy_approximation) == 0) {
      NA_real_
    } else {
      cor(comparison$occupancy_simulation,
          comparison$occupancy_approximation)
    },
    RMSE_first_500 = sqrt(mean(early$residual^2)),
    max_absolute_deviation_first_500 = max(abs(early$residual)),
    simulation_time_to_90pct = first_time_at_or_above(
      comparison$time_post_burnout, comparison$occupancy_simulation, 0.9
    ),
    approximation_time_to_90pct = first_time_at_or_above(
      comparison$time_post_burnout, comparison$occupancy_approximation, 0.9
    ),
    simulation_time_to_95pct = first_time_at_or_above(
      comparison$time_post_burnout, comparison$occupancy_simulation, 0.95
    ),
    approximation_time_to_95pct = first_time_at_or_above(
      comparison$time_post_burnout, comparison$occupancy_approximation, 0.95
    ),
    stringsAsFactors = FALSE
  )

  list(curves = curves, diagnostics = diagnostics, params = params)
}

analyses <- lapply(seq_along(metapop_results), function(i) {
  analyse_one(metapop_results[[i]], established_results[[i]], i)
})
all_curves <- bind_rows(lapply(analyses, `[[`, "curves"))
all_diagnostics <- bind_rows(lapply(analyses, `[[`, "diagnostics"))

palette <- c(
  "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"
)
groups <- unique(all_diagnostics$group)

for (group_name in groups) {
  group_rows <- all_diagnostics |>
    filter(group == group_name) |>
    arrange(result_id)
  label_levels <- vapply(
    seq_len(nrow(group_rows)),
    function(i) format_parameter(
      group_name, group_rows$parameter_value[i]
    ),
    character(1)
  )
  plot_data <- all_curves |>
    filter(group == group_name) |>
    mutate(parameter_label = factor(parameter_label, levels = label_levels))

  fixed_params <- analyses[[group_rows$result_id[1]]]$params
  bounded_note <- if (any(group_rows$P1_was_bounded)) {
    " P1 estimates at 0 or 1 were bounded for numerical calculation."
  } else {
    ""
  }

  p <- ggplot(
    plot_data,
    aes(
      time_post_burnout,
      occupancy,
      colour = parameter_label,
      linetype = trajectory,
      group = interaction(parameter_label, trajectory)
    )
  ) +
    geom_line(linewidth = 0.65) +
    scale_colour_manual(values = palette[seq_along(label_levels)]) +
    scale_linetype_manual(values = c(
      "Established stochastic occupancy (tau = 50)" = "solid",
      "Analytical approximation" = "22"
    )) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      labels = scales::percent_format(accuracy = 1)
    ) +
    labs(
      x = "Time since each raw-occupancy minimum",
      y = "Fraction of established patches (tau = 50)",
      colour = sprintf("Varied parameter: %s", group_name),
      linetype = NULL,
      title = sprintf(
        "Patch occupancy approximation: varying %s", group_name
      ),
      subtitle = fixed_parameter_subtitle(group_name, fixed_params),
      caption = paste0(
        "Solid curves use tau=50 established stochastic occupancy.",
        " Analytical time zero is the corresponding raw-occupancy minimum.",
        bounded_note
      )
    ) +
    theme(
      legend.position = "right",
      plot.caption = element_text(hjust = 0, size = 8)
    )

  plot_file <- file.path(
    outdir,
    sprintf("one_state_established_tau50_compare_%s.pdf", group_name)
  )
  tryCatch(
    ggsave(plot_file, p, width = 11, height = 6.5),
    error = function(e) {
      warning(
        "Could not write ", plot_file,
        " (it may be open in another program): ", conditionMessage(e)
      )
    }
  )
}

write.csv(all_diagnostics, diagnostic_file, row.names = FALSE)

cat("Single-patch P1 source: R0 x K x r extinction grid at t=200; ",
    "Poisson initial infected count with mean 10.\n", sep = "")
cat("Exact grid values are preferred; non-grid R0/K values use a matching-r ",
    "two-dimensional GAM.\n", sep = "")
cat("Stochastic target: established occupancy with tau=50; alignment: raw occupancy minimum.\n")
cat("Probability bound used for analytical calculation: ",
    probability_tolerance, " to ", 1 - probability_tolerance, "\n", sep = "")
cat("Scenarios analysed: ", nrow(all_diagnostics), "\n", sep = "")
cat("P1 values bounded for numerical calculation: ",
    sum(all_diagnostics$P1_was_bounded), "\n", sep = "")
cat("Median RMSE: ", median(all_diagnostics$RMSE), "\n", sep = "")
cat("Median first-500 RMSE: ", median(all_diagnostics$RMSE_first_500), "\n", sep = "")
cat("Outputs: ", outdir, "\n", sep = "")
