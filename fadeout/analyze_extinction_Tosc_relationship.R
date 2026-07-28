## Relationship between deterministic oscillation period and the conditional
## mean stochastic extinction time in the existing one-patch demography grid.
##
## Run from the repository root:
##   Rscript fadeout/analyze_extinction_Tosc_relationship.R
##
## This script only reads the combined SHARCNET output. It does not run any
## stochastic simulations.

library(dplyr)
library(ggplot2)
library(here)

theme_set(theme_bw())

input_file <- here::here(
  "odin", "sharcnet", "outputs",
  "euler_onepatch_onestrain_extinct_logistic_continuous_demoggrid.rds"
)
run_script <- here::here(
  "odin", "sharcnet", "euler_onepatch_onestrain_extinct_run_array.R"
)
summary_script <- here::here("plagueMetapop", "R", "discrete_run.R")
two_state_script <- here::here("fadeout", "two_state_functions.R")

outdir <- here::here("fadeout", "output", "extinction_Tosc_relationship")
datadir <- file.path(outdir, "data")
figdir <- file.path(outdir, "figures")
dir.create(datadir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

grid_file <- file.path(datadir, "extinction_Tosc_grid.csv")
regression_file <- file.path(datadir, "extinction_Tosc_regression.csv")
ratio_file <- file.path(datadir, "extinction_Tosc_ratio_summary.csv")
main_figure <- file.path(figdir, "extinction_time_vs_Tosc_loglog.pdf")
ratio_figure <- file.path(figdir, "extinction_Tosc_ratio.pdf")

required_files <- c(input_file, run_script, summary_script, two_state_script)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Required input/source files not found: ",
       paste(missing_files, collapse = ", "))
}

## Metadata verified below against the actual SHARCNET array script.
gamma <- 1
dt <- 0.1
nsim <- 1000L
t_max <- 200
T_osc_max <- 50
n_patch <- 1L
alpha <- 0
I_init <- 10

grid <- readRDS(input_file)
required_columns <- c(
  "R0", "K", "r", "mean_ext_time.I1", "ext_prob.I1"
)
if (!is.data.frame(grid)) stop("Input RDS is not a data frame")
missing_columns <- setdiff(required_columns, names(grid))
if (length(missing_columns)) {
  stop("Input grid is missing columns: ", paste(missing_columns, collapse = ", "))
}

## Validate the metadata against executable source rather than relying only on
## the README. Whitespace is removed to make these literal checks robust.
run_text <- paste(readLines(run_script, warn = FALSE), collapse = "")
run_compact <- gsub("[[:space:]]+", "", run_text)
source_checks <- c(
  nsim = grepl("nsim<-1000L", run_compact, fixed = TRUE),
  dt = grepl("dt<-if(reedfrost==1)1else0.1", run_compact, fixed = TRUE),
  t_max = grepl("nt=round(200/dt)", run_compact, fixed = TRUE),
  n_patch = grepl("n_patch=1", run_compact, fixed = TRUE),
  alpha = grepl("alpha=0", run_compact, fixed = TRUE),
  I_init = grepl("I_init=c(10,0)", run_compact, fixed = TRUE),
  gamma = grepl("gamma=c(1,1)", run_compact, fixed = TRUE),
  logistic = grepl("logistic_growth=logistic_growth", run_compact, fixed = TRUE),
  continuous = grepl("dt<-if(reedfrost==1)1else0.1", run_compact, fixed = TRUE)
)
if (!all(source_checks)) {
  stop("Simulation metadata validation failed for: ",
       paste(names(source_checks)[!source_checks], collapse = ", "))
}

## sumfun_discrete() defines ext1 as the row index of the first saved row with
## total I1 == 0. run_simulator_odin() retains the initial-condition row at
## step 0, and successive saved rows are dt apart. Therefore row index j maps
## to model time (j - 1) * dt. The stored mean is conditional on non-NA ext1,
## i.e. conditional on extinction within the finite observation window.
summary_text <- paste(readLines(summary_script, warn = FALSE), collapse = "")
summary_compact <- gsub("[[:space:]]+", "", summary_text)
index_conversion_verified <-
  grepl("first_ext<-function", summary_compact, fixed = TRUE) &&
  grepl("idx[1L]", summary_compact, fixed = TRUE) &&
  grepl("mean_ext_time.I1=mean(ext1,na.rm=TRUE)",
        summary_compact, fixed = TRUE)
if (!index_conversion_verified) {
  stop("Could not verify mean_ext_time.I1 row-index semantics")
}

expected_R0 <- seq(1.1, 5, by = 0.1)
expected_K <- 10^seq(3, 6, by = 0.25)
expected_r <- c(0.05, 0.1, 0.125, 0.2, 0.4)

same_numeric_set <- function(observed, expected, tolerance = 1e-10) {
  observed <- sort(unique(observed))
  expected <- sort(expected)
  length(observed) == length(expected) &&
    all(abs(observed - expected) <=
          tolerance * pmax(1, abs(expected)))
}

grid_checks <- c(
  rows_2600 = nrow(grid) == 2600L,
  R0_grid = same_numeric_set(grid$R0, expected_R0),
  K_grid = same_numeric_set(grid$K, expected_K),
  r_grid = same_numeric_set(grid$r, expected_r),
  no_duplicate_cells = !anyDuplicated(grid[c("R0", "K", "r")])
)
if (!all(grid_checks)) {
  stop("Grid validation failed for: ",
       paste(names(grid_checks)[!grid_checks], collapse = ", "))
}

omega_sq <- gamma * grid$r * (grid$R0 - 1) -
  grid$r^2 * grid$R0^2 / 4
oscillatory <- is.finite(omega_sq) & omega_sq > 0
T_osc <- rep(NA_real_, nrow(grid))
T_osc[oscillatory] <- 2 * pi / sqrt(omega_sq[oscillatory])

## Secondary sensitivity diagnostic from the exact local Jacobian of the
## deterministic logistic model. It is retained in the table but is not used
## in the main regression or figures.
omega_logistic_sq <- gamma * grid$r * (grid$R0 - 1) / grid$R0 -
  grid$r^2 / (4 * grid$R0^2)
T_osc_logistic <- rep(NA_real_, nrow(grid))
logistic_oscillatory <- is.finite(omega_logistic_sq) &
  omega_logistic_sq > 0
T_osc_logistic[logistic_oscillatory] <-
  2 * pi / sqrt(omega_logistic_sq[logistic_oscillatory])

mean_T_ext <- (grid$mean_ext_time.I1 - 1) * dt
positive_extinction_time <- is.finite(mean_T_ext) & mean_T_ext > 0
within_Tosc_window <- is.finite(T_osc) & T_osc > 0 & T_osc < T_osc_max
finite_valid <- oscillatory & within_Tosc_window & positive_extinction_time

analysis_grid <- grid |>
  transmute(
    R0, K, log10K = log10(K), r,
    gamma = gamma, dt = dt, nsim = nsim, t_max = t_max,
    n_patch = n_patch, alpha = alpha, I_init = I_init,
    ext_prob = ext_prob.I1,
    n_ext = round(nsim * ext_prob.I1),
    mean_ext_time_index = mean_ext_time.I1,
    mean_T_ext = mean_T_ext,
    omega_sq = omega_sq,
    T_osc = T_osc,
    omega_logistic_sq = omega_logistic_sq,
    T_osc_logistic = T_osc_logistic,
    ratio = if_else(finite_valid, mean_T_ext / T_osc, NA_real_),
    log10_T_osc = if_else(finite_valid, log10(T_osc), NA_real_),
    log10_mean_T_ext =
      if_else(finite_valid, log10(mean_T_ext), NA_real_),
    log10_ratio =
      if_else(finite_valid, log10(mean_T_ext / T_osc), NA_real_),
    oscillatory = oscillatory,
    logistic_oscillatory = logistic_oscillatory,
    positive_extinction_time = positive_extinction_time,
    within_Tosc_window = within_Tosc_window,
    finite_valid = finite_valid
  )

valid <- analysis_grid |> filter(finite_valid)
if (nrow(valid) < 3L) stop("Too few finite valid grid cells for regression")

fit <- lm(log10_mean_T_ext ~ log10_T_osc, data = valid)
fit_summary <- summary(fit)
fit_ci <- confint(fit, level = 0.95)
coef_table <- as.data.frame(fit_summary$coefficients)
coef_table$term <- rownames(coef_table)
rownames(coef_table) <- NULL
names(coef_table)[1:4] <- c("estimate", "std_error", "t_value", "p_value")
coef_table$conf_low <- fit_ci[coef_table$term, 1]
coef_table$conf_high <- fit_ci[coef_table$term, 2]

## Quantify, but do not model, the vertical K spread at a fixed T_osc
## (equivalently at fixed R0 and r).
within_column_spread <- valid |>
  group_by(R0, r, T_osc) |>
  summarise(
    n_K = n(),
    min_mean_T_ext = min(mean_T_ext),
    max_mean_T_ext = max(mean_T_ext),
    fold_spread = max_mean_T_ext / min_mean_T_ext,
    log10_range = log10(max_mean_T_ext) - log10(min_mean_T_ext),
    .groups = "drop"
  ) |>
  filter(n_K >= 2L)

regression_output <- coef_table |>
  mutate(
    r_squared = fit_summary$r.squared,
    adjusted_r_squared = fit_summary$adj.r.squared,
    residual_standard_error = fit_summary$sigma,
    n_grid_cells = nrow(valid),
    T_osc_max_exclusive = T_osc_max,
    implied_c_from_intercept = 10^coef(fit)[["(Intercept)"]],
    median_within_column_fold_spread =
      median(within_column_spread$fold_spread),
    p95_within_column_fold_spread =
      unname(quantile(within_column_spread$fold_spread, 0.95)),
    maximum_within_column_fold_spread =
      max(within_column_spread$fold_spread)
  ) |>
  select(
    term, estimate, std_error, conf_low, conf_high, t_value, p_value,
    r_squared, adjusted_r_squared, residual_standard_error, n_grid_cells,
    T_osc_max_exclusive, implied_c_from_intercept,
    median_within_column_fold_spread,
    p95_within_column_fold_spread, maximum_within_column_fold_spread
  )

ratio_values <- valid$ratio
ratio_summary <- data.frame(
  n_grid_cells = length(ratio_values),
  T_osc_max_exclusive = T_osc_max,
  minimum = min(ratio_values),
  percentile_05 = unname(quantile(ratio_values, 0.05)),
  percentile_25 = unname(quantile(ratio_values, 0.25)),
  median = median(ratio_values),
  geometric_mean = exp(mean(log(ratio_values))),
  percentile_75 = unname(quantile(ratio_values, 0.75)),
  percentile_95 = unname(quantile(ratio_values, 0.95)),
  maximum = max(ratio_values),
  median_within_column_fold_spread =
    median(within_column_spread$fold_spread),
  percentile_95_within_column_fold_spread =
    unname(quantile(within_column_spread$fold_spread, 0.95)),
  maximum_within_column_fold_spread =
    max(within_column_spread$fold_spread)
)

## Cross-check the baseline formula against the actual two-state helper.
source(two_state_script)
baseline_summary <- compute_transient_outbreak_summary(
  R0 = 2.5, K = 10000, r = 0.125, gamma = gamma,
  I0 = 10, c = 0.5
)
baseline_row <- analysis_grid |>
  filter(abs(R0 - 2.5) < 1e-10, abs(K - 10000) < 1e-8,
         abs(r - 0.125) < 1e-10)
baseline_present <- nrow(baseline_row) == 1L
formula_agrees <- baseline_present &&
  isTRUE(all.equal(
    baseline_row$omega_sq, baseline_summary$exact_argument,
    tolerance = 1e-12
  )) &&
  isTRUE(all.equal(
    baseline_row$T_osc, baseline_summary$T_osc,
    tolerance = 1e-12
  ))

## Explicit numerical validations.
period_valid <- all(
  is.na(analysis_grid$T_osc[!analysis_grid$oscillatory])
) && all(
  is.finite(analysis_grid$T_osc[analysis_grid$oscillatory]) &
    analysis_grid$T_osc[analysis_grid$oscillatory] > 0
)
derived_valid <- all(
  is.na(analysis_grid$ratio[!analysis_grid$finite_valid])
) && all(
  is.finite(analysis_grid$ratio[analysis_grid$finite_valid]) &
    analysis_grid$ratio[analysis_grid$finite_valid] > 0 &
    is.finite(analysis_grid$log10_T_osc[analysis_grid$finite_valid]) &
    is.finite(analysis_grid$log10_mean_T_ext[
      analysis_grid$finite_valid
    ])
)
same_period_across_K <- analysis_grid |>
  group_by(R0, r) |>
  summarise(n_period = n_distinct(T_osc, na.rm = TRUE), .groups = "drop") |>
  summarise(ok = all(n_period <= 1L)) |>
  pull(ok)

validation <- c(
  grid_checks,
  metadata_matches_source = all(source_checks),
  extinction_index_conversion_verified = index_conversion_verified,
  period_positive_only_when_oscillatory = period_valid,
  ratios_and_logs_only_for_valid_rows = derived_valid,
  same_Tosc_across_K = same_period_across_K,
  baseline_present = baseline_present,
  formula_matches_two_state_helper = formula_agrees
)
if (!all(validation)) {
  stop("Analysis validation failed for: ",
       paste(names(validation)[!validation], collapse = ", "))
}

write.csv(analysis_grid, grid_file, row.names = FALSE)
write.csv(regression_output, regression_file, row.names = FALSE)
write.csv(ratio_summary, ratio_file, row.names = FALSE)

facet_levels <- as.character(expected_r)
plot_data <- valid |>
  mutate(r_facet = factor(as.character(r), levels = facet_levels))

x_range <- range(valid$T_osc)
x_values <- 10^seq(log10(x_range[1]), log10(x_range[2]), length.out = 200)
reference_data <- expand.grid(
  T_osc = x_values,
  c_value = c(0.25, 0.5, 1, 2),
  r_facet = factor(facet_levels, levels = facet_levels)
) |>
  mutate(
    mean_T_ext = c_value * T_osc,
    reference = factor(
      paste0("c = ", c_value),
      levels = c("c = 0.25", "c = 0.5", "c = 1", "c = 2")
    )
  )
regression_data <- expand.grid(
  T_osc = x_values,
  r_facet = factor(facet_levels, levels = facet_levels)
) |>
  mutate(
    mean_T_ext = 10^predict(
      fit, newdata = data.frame(log10_T_osc = log10(T_osc))
    )
  )

reference_colours <- c(
  "c = 0.25" = "grey65", "c = 0.5" = "#D55E00",
  "c = 1" = "grey40", "c = 2" = "grey20"
)
reference_widths <- c(
  "c = 0.25" = 0.45, "c = 0.5" = 0.9, "c = 1" = 0.45, "c = 2" = 0.45
)

p_main <- ggplot(plot_data, aes(T_osc, mean_T_ext)) +
  geom_line(
    data = reference_data,
    aes(colour = reference, linewidth = reference),
    alpha = 0.9
  ) +
  geom_line(
    data = regression_data,
    colour = "black", linewidth = 0.9, linetype = "dashed"
  ) +
  geom_point(aes(fill = log10K), shape = 21, size = 1.7, alpha = 0.62,
             colour = "black", stroke = 0.12) +
  facet_wrap(
    ~r_facet,
    labeller = labeller(r_facet = function(x) paste0("r = ", x))
  ) +
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_viridis_c(option = "C", name = expression(log[10](K))) +
  scale_colour_manual(values = reference_colours, name = "Proportional lines") +
  scale_linewidth_manual(values = reference_widths, guide = "none") +
  labs(
    x = expression("Deterministic oscillation period " * T[osc]),
    y = "Conditional mean extinction time among extinct runs",
    title = "Conditional mean extinction time versus deterministic oscillation period",
    subtitle = paste0(
      "Each point is one R0-K-r grid cell; K is not included in the regression\n",
      "Dashed black line: overall unweighted log-log regression; T_osc < 50"
    ),
    caption = paste0(
      "T_osc depends on R0 and r under the current two-state formula, not K; ",
      "different K values therefore share x-coordinates.\n",
      "All grid cells are separate regression observations. Response is ",
      "conditional on extinction within t_max = 200; only T_osc < 50 is analysed."
    )
  ) +
  guides(
    fill = guide_colourbar(order = 1),
    colour = guide_legend(order = 2)
  ) +
  theme(
    legend.position = "right",
    plot.caption = element_text(hjust = 0, size = 8),
    plot.subtitle = element_text(size = 9)
  )

ggsave(main_figure, p_main, width = 12, height = 8)

p_ratio <- ggplot(plot_data, aes(T_osc, ratio)) +
  geom_hline(yintercept = c(0.25, 1, 2), colour = "grey45",
             linewidth = 0.45, linetype = "dashed") +
  geom_hline(yintercept = 0.5, colour = "#D55E00", linewidth = 0.9) +
  geom_point(aes(fill = log10K), shape = 21, size = 1.8, alpha = 0.65,
             colour = "black", stroke = 0.12) +
  facet_wrap(
    ~r_facet,
    labeller = labeller(r_facet = function(x) paste0("r = ", x))
  ) +
  scale_x_log10() +
  scale_y_log10() +
  scale_fill_viridis_c(option = "C", name = expression(log[10](K))) +
  labs(
    x = expression("Deterministic oscillation period " * T[osc]),
    y = expression("Conditional mean extinction-time ratio " *
                     mean(T[ext]) / T[osc]),
    title = "Conditional extinction time relative to oscillation period",
    subtitle = paste0(
      "Each point is one R0-K-r grid cell; orange line marks current c = 0.5; ",
      "T_osc < 50"
    ),
    caption = paste0(
      "Grey dashed reference lines: c = 0.25, 1, and 2. ",
      "K is shown descriptively and is not fitted."
    )
  ) +
  theme(
    legend.position = "right",
    plot.caption = element_text(hjust = 0, size = 8),
    plot.subtitle = element_text(size = 9)
  )

ggsave(ratio_figure, p_ratio, width = 12, height = 8)

cat("Validation checks passed: ", length(validation), "/", length(validation),
    "\n", sep = "")
cat("Grid rows: ", nrow(analysis_grid), "; regression rows: ", nrow(valid),
    "; non-oscillatory rows: ", sum(!analysis_grid$oscillatory),
    "; rows with T_osc >= ", T_osc_max, ": ",
    sum(analysis_grid$oscillatory & !analysis_grid$within_Tosc_window),
    "; invalid/non-positive extinction-time rows: ",
    sum(!analysis_grid$positive_extinction_time), "\n", sep = "")
cat("Extinction-time conversion: (mean_ext_time.I1 - 1) * dt, dt = ",
    dt, "\n", sep = "")
cat("Slope: ", coef(fit)[["log10_T_osc"]], " (95% CI ",
    fit_ci["log10_T_osc", 1], ", ", fit_ci["log10_T_osc", 2], ")\n",
    sep = "")
cat("Intercept: ", coef(fit)[["(Intercept)"]], "; 10^intercept = ",
    10^coef(fit)[["(Intercept)"]], "; R-squared = ",
    fit_summary$r.squared, "\n", sep = "")
cat("Ratio median: ", ratio_summary$median, "; 25%-75%: ",
    ratio_summary$percentile_25, "-", ratio_summary$percentile_75,
    "; 5%-95%: ", ratio_summary$percentile_05, "-",
    ratio_summary$percentile_95, "\n", sep = "")
cat("Within-(R0,r) K fold spread: median ",
    ratio_summary$median_within_column_fold_spread, "; 95th percentile ",
    ratio_summary$percentile_95_within_column_fold_spread, "; maximum ",
    ratio_summary$maximum_within_column_fold_spread, "\n", sep = "")
cat("Outputs saved under: ", outdir, "\n", sep = "")
