## Small local validation suite for the logistic burnout approximation.
## Run from the repository root:
##   Rscript fadeout/logistic_burnout/validate_logistic_burnout.R

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

## A. Closed-form logistic recovery versus numerical ODE integration.
recovery_cases <- expand.grid(
  x_in = c(0.02, 0.2, 0.7, 0.99),
  r = c(0.02, 0.125, 0.3)
)
recovery_results <- do.call(rbind, lapply(
  seq_len(nrow(recovery_cases)),
  function(i) {
    p <- recovery_cases[i, ]
    times <- seq(0, 100, length.out = 301)
    numerical <- as.data.frame(deSolve::ode(
      y = c(x = p$x_in), times = times,
      func = function(t, state, parameters) {
        list(c(x = p$r * state[["x"]] * (1 - state[["x"]])))
      },
      parms = NULL, rtol = 1e-11, atol = 1e-13
    ))
    closed <- logistic_recovery_x(times, p$x_in, p$r)
    data.frame(
      test = "logistic_recovery", R0 = NA_real_, r = p$r,
      x_in = p$x_in, time = NA_real_, analytical = NA_real_,
      numerical = NA_real_,
      absolute_difference = max(abs(closed - numerical$x)),
      tolerance = 2e-8, passed = max(abs(closed - numerical$x)) < 2e-8,
      integration_converged = NA, integration_message = NA_character_,
      q1 = NA_real_, q1_tight = NA_real_,
      simulated_q1 = NA_real_, monte_carlo_se = NA_real_,
      n_replicates = NA_integer_, n_time_censored = NA_integer_
    )
  }
))

## B. Analytical H(t) versus direct numerical quadrature.
H_cases <- expand.grid(
  R0 = c(1.2, 2.5, 5),
  r = c(0.02, 0.125, 0.3),
  x_in = c(0.05, 0.3),
  time = c(1, 10, 50)
)
H_results <- do.call(rbind, lapply(
  seq_len(nrow(H_cases)),
  function(i) {
    p <- H_cases[i, ]
    analytical <- logistic_lineage_H(p$time, p$R0, p$r, p$x_in)
    numerical <- stats::integrate(
      function(t) p$R0 * logistic_recovery_x(t, p$x_in, p$r) - 1,
      0, p$time, rel.tol = 1e-11
    )$value
    difference <- abs(analytical - numerical)
    data.frame(
      test = "integrated_net_growth", R0 = p$R0, r = p$r,
      x_in = p$x_in, time = p$time, analytical = analytical,
      numerical = numerical, absolute_difference = difference,
      tolerance = 2e-8, passed = difference < 2e-8,
      integration_converged = NA, integration_message = NA_character_,
      q1 = NA_real_, q1_tight = NA_real_,
      simulated_q1 = NA_real_, monte_carlo_se = NA_real_,
      n_replicates = NA_integer_, n_time_censored = NA_integer_
    )
  }
))

## C. Improper-integral convergence at deterministic boundary entries.
integral_cases <- expand.grid(
  R0 = c(1.2, 2, 3, 5),
  r = c(0.02, 0.05, 0.125, 0.3)
)
integral_results <- do.call(rbind, lapply(
  seq_len(nrow(integral_cases)),
  function(i) {
    p <- integral_cases[i, ]
    burnout <- logistic_burnout_probability(
      p$R0, p$r, K = 10000, I0 = 1,
      rel.tol = 1e-8, subdivisions = 500
    )
    tight <- if (burnout$boundary_entry_found) {
      logistic_lineage_extinction_probability(
        p$R0, p$r, burnout$x_in,
        rel.tol = 1e-10, subdivisions = 1000
      )
    } else {
      list(q1 = NA_real_, converged = FALSE, message = burnout$status)
    }
    difference <- abs(burnout$q1 - tight$q1)
    passed <- burnout$boundary_entry_found &&
      burnout$integration_converged && tight$converged &&
      burnout$q1 >= 0 && burnout$q1 <= 1 &&
      difference < 2e-7
    data.frame(
      test = "extinction_integral", R0 = p$R0, r = p$r,
      x_in = burnout$x_in, time = NA_real_, analytical = NA_real_,
      numerical = NA_real_, absolute_difference = difference,
      tolerance = 2e-7, passed = passed,
      integration_converged = burnout$integration_converged &&
        tight$converged,
      integration_message = paste(
        burnout$integration_message, tight$message, sep = " | "
      ),
      q1 = burnout$q1, q1_tight = tight$q1,
      simulated_q1 = NA_real_, monte_carlo_se = NA_real_,
      n_replicates = NA_integer_, n_time_censored = NA_integer_
    )
  }
))

## D. Exact thinning validation of the nonhomogeneous branching process.
branching_cases <- data.frame(
  R0 = c(1.2, 2, 3, 5),
  r = c(0.05, 0.125, 0.125, 0.3)
)
branching_results <- do.call(rbind, lapply(
  seq_len(nrow(branching_cases)),
  function(i) {
    p <- branching_cases[i, ]
    burnout <- logistic_burnout_probability(
      p$R0, p$r, K = 10000, I0 = 1
    )
    if (!burnout$boundary_entry_found || !burnout$integration_converged) {
      return(data.frame(
        test = "branching_simulation", R0 = p$R0, r = p$r,
        x_in = burnout$x_in, time = NA_real_, analytical = NA_real_,
        numerical = NA_real_, absolute_difference = NA_real_,
        tolerance = NA_real_, passed = FALSE,
        integration_converged = burnout$integration_converged,
        integration_message = burnout$status, q1 = burnout$q1,
        q1_tight = NA_real_, simulated_q1 = NA_real_,
        monte_carlo_se = NA_real_, n_replicates = NA_integer_,
        n_time_censored = NA_integer_
      ))
    }
    simulation <- simulate_logistic_lineage(
      p$R0, p$r, burnout$x_in, n_replicates = 2500,
      seed = 20260730 + i, tmax = 500, population_cap = 300
    )
    difference <- abs(simulation$extinction_frequency - burnout$q1)
    tolerance <- 4 * simulation$monte_carlo_se + 0.01
    data.frame(
      test = "branching_simulation", R0 = p$R0, r = p$r,
      x_in = burnout$x_in, time = NA_real_, analytical = burnout$q1,
      numerical = simulation$extinction_frequency,
      absolute_difference = difference, tolerance = tolerance,
      passed = difference <= tolerance &&
        simulation$n_time_censored == 0,
      integration_converged = burnout$integration_converged,
      integration_message = burnout$integration_message,
      q1 = burnout$q1, q1_tight = NA_real_,
      simulated_q1 = simulation$extinction_frequency,
      monte_carlo_se = simulation$monte_carlo_se,
      n_replicates = simulation$n_replicates,
      n_time_censored = simulation$n_time_censored
    )
  }
))

validation <- rbind(
  recovery_results, H_results, integral_results, branching_results
)
write.csv(
  validation, file.path(output_dir, "validation_results.csv"),
  row.names = FALSE
)

plot_data <- branching_results
p <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = q1, y = simulated_q1,
    colour = factor(R0), shape = factor(r)
  )
) +
  ggplot2::geom_abline(
    slope = 1, intercept = 0, linetype = "dashed", colour = "grey40"
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = pmax(0, simulated_q1 - 1.96 * monte_carlo_se),
      ymax = pmin(1, simulated_q1 + 1.96 * monte_carlo_se)
    ),
    width = 0
  ) +
  ggplot2::geom_point(size = 3) +
  ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  ggplot2::labs(
    x = "Analytical-numerical single-lineage extinction probability q1",
    y = "Thinning-simulation extinction frequency",
    colour = "R0", shape = "r",
    title = "Validation of the logistic boundary-layer branching process",
    subtitle = paste0(
      "2,500 replicates per case; exact thinning envelope; ",
      "K = 10,000 and I(0) = 1 used only to obtain x_in"
    ),
    caption = paste0(
      "Bars are approximate 95% Monte Carlo intervals. Surviving lineages ",
      "are followed to I = 300 or 500 disease generations."
    )
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    legend.position = "right",
    plot.subtitle = ggplot2::element_text(size = 9),
    plot.caption = ggplot2::element_text(hjust = 0, size = 8)
  )
ggplot2::ggsave(
  file.path(figure_dir, "validation_q1.png"),
  p, width = 8, height = 7, dpi = 180
)

if (!all(validation$passed)) {
  failed <- validation[!validation$passed, ]
  stop(
    "Validation failures: ",
    paste(unique(failed$test), collapse = ", "),
    ". Inspect outputs/validation_results.csv"
  )
}
cat("Validation checks passed: ", nrow(validation), "/", nrow(validation),
    "\n", sep = "")
print(branching_results[c(
  "R0", "r", "x_in", "q1", "simulated_q1",
  "monte_carlo_se", "n_time_censored"
)])
