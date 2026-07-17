## Compare the analytical patch-occupancy approximation with every existing
## trajectory in the one-factor-at-a-time stochastic metapopulation analysis.
## This script reads existing outputs only; it does not run simulations.

library(dplyr)
library(ggplot2)
library(here)
library(mgcv)
library(plagueMetapop)

theme_set(theme_bw())

single_patch_file <- here::here(
  "odin", "sharcnet", "outputs",
  "euler_onepatch_onestrain_extinct_logistic_continuous.rds"
)
metapop_file <- here::here(
  "fadeout", "output", "patch_occupancy", "data",
  "occupancy_results.rds"
)
outdir <- here::here("fadeout", "output", "patch_occupancy_approximation")
diagnostic_file <- file.path(outdir, "patch_occupancy_approximation_diagnostics.csv")

for (filename in c(single_patch_file, metapop_file)) {
  if (!file.exists(filename)) stop("Required input file not found: ", filename)
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

single_patch <- readRDS(single_patch_file)
required_single <- c("R0", "K", "ext_prob.I1")
missing_single <- setdiff(required_single, names(single_patch))
if (!is.data.frame(single_patch) || nrow(single_patch) == 0L) {
  stop("Single-patch input must be a non-empty data frame")
}
if (length(missing_single) > 0L) {
  stop("Single-patch data are missing columns: ",
       paste(missing_single, collapse = ", "))
}

metapop_results <- readRDS(metapop_file)
if (!is.list(metapop_results) || length(metapop_results) == 0L) {
  stop("Metapopulation input must be a non-empty list")
}

## Follow notes/notes_16jul.rmd. This emulator is calibrated at the r value
## used in the single-patch grid (r=0.125); applying it across the r comparison
## assumes P1 depends only on R0 and K, as in the proposed approximation.
gam_fit <- mgcv::gam(
  ext_prob.I1 ~ te(R0, K, k = c(12, 12)),
  data = single_patch,
  method = "REML"
)

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
  paste(vapply(
    fixed,
    function(name) format_parameter(name, params[[name]]),
    character(1)
  ), collapse = "; ")
}

first_time_at_or_above <- function(time, value, level) {
  index <- which(value >= level)[1L]
  if (is.na(index)) NA_real_ else time[index]
}

analyse_one <- function(result, result_id) {
  required_result <- c(
    "varied_parameter", "parameter_value", "params", "raw_meta_summary"
  )
  missing_result <- setdiff(required_result, names(result))
  if (length(missing_result) > 0L) {
    stop("Metapopulation result ", result_id, " is missing fields: ",
         paste(missing_result, collapse = ", "))
  }

  params <- result$params
  simulation <- result$raw_meta_summary
  required_meta <- c("step", "occupied_patches", "occupancy_fraction")
  missing_meta <- setdiff(required_meta, names(simulation))
  if (length(missing_meta) > 0L) {
    stop("Trajectory ", result_id, " is missing columns: ",
         paste(missing_meta, collapse = ", "))
  }

  extinction_prediction <- as.numeric(predict(
    gam_fit,
    newdata = data.frame(R0 = params$R0, K = params$K),
    type = "response"
  ))
  P1_GAM_raw <- 1 - extinction_prediction
  P1 <- pmin(pmax(P1_GAM_raw, probability_tolerance),
             1 - probability_tolerance)
  P1_was_bounded <- !isTRUE(all.equal(P1_GAM_raw, P1))

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
    stop("Invalid endemic I* for result ", result_id, ": ", I_star)
  }

  ## The analytical initial epidemic is instantaneous. Align t=0 with the
  ## minimum raw simulated occupancy. If the trajectory never falls below one,
  ## which.min() returns simulation time zero and the diagnostic records this.
  minimum_index <- which.min(simulation$occupancy_fraction)
  time_shift <- simulation$step[minimum_index]
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
        trajectory = "Existing stochastic trajectory"
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
    P1_GAM_raw = P1_GAM_raw,
    P1_used = P1,
    P1_was_bounded = P1_was_bounded,
    P1_exact_grid = P1_exact_grid,
    I_star = I_star,
    lambda = lambda,
    simulation_minimum = simulation$occupancy_fraction[minimum_index],
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
  analyse_one(metapop_results[[i]], i)
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
    " GAM predictions outside [0,1] were bounded for calculation."
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
      "Existing stochastic trajectory" = "solid",
      "Analytical approximation" = "22"
    )) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      labels = scales::percent_format(accuracy = 1)
    ) +
    labs(
      x = "Time since each simulated occupancy minimum",
      y = "Fraction of occupied patches",
      colour = sprintf("Varied parameter: %s", group_name),
      linetype = NULL,
      title = sprintf(
        "Patch occupancy approximation: varying %s", group_name
      ),
      subtitle = fixed_parameter_subtitle(group_name, fixed_params),
      caption = paste0(
        "Each curve uses one existing stochastic trajectory; no uncertainty interval.",
        " Analytical time zero is the corresponding simulated occupancy minimum.",
        bounded_note
      )
    ) +
    theme(
      legend.position = "right",
      plot.caption = element_text(hjust = 0, size = 8)
    )

  plot_file <- file.path(
    outdir,
    sprintf("patch_occupancy_approximation_compare_%s.pdf", group_name)
  )
  ggsave(plot_file, p, width = 11, height = 6.5)
}

write.csv(all_diagnostics, diagnostic_file, row.names = FALSE)

cat("Single-patch P1 source: GAM fitted to extinction by t=200; ",
    "Poisson initial infected count with mean 10.\n", sep = "")
cat("P1 emulator calibration uses r=0.125; the r comparison assumes P1(R0,K).\n")
cat("Probability bound used for analytical calculation: ",
    probability_tolerance, " to ", 1 - probability_tolerance, "\n", sep = "")
cat("Scenarios analysed: ", nrow(all_diagnostics), "\n", sep = "")
cat("GAM predictions bounded: ", sum(all_diagnostics$P1_was_bounded), "\n", sep = "")
cat("Median RMSE: ", median(all_diagnostics$RMSE), "\n", sep = "")
cat("Median first-500 RMSE: ", median(all_diagnostics$RMSE_first_500), "\n", sep = "")
cat("Outputs: ", outdir, "\n", sep = "")
