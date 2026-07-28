## One-parameter-at-a-time visualization of the existing one-patch,
## one-strain stochastic extinction grid.
##
## Run from the repository root:
##   Rscript fadeout/analyze_onepatch_extinction_sensitivity.R
##
## This script reads the combined SHARCNET summary and never runs simulations.

library(dplyr)
library(ggplot2)
library(here)
library(cowplot)

theme_set(theme_bw())

input_file <- here::here(
  "odin", "sharcnet", "outputs",
  "euler_onepatch_onestrain_extinct_logistic_continuous_demoggrid.rds"
)
run_script <- here::here(
  "odin", "sharcnet", "euler_onepatch_onestrain_extinct_run_array.R"
)
summary_script <- here::here("plagueMetapop", "R", "discrete_run.R")

outdir <- here::here("fadeout", "output", "onepatch_extinction_sensitivity")
datadir <- file.path(outdir, "data")
figdir <- file.path(outdir, "figures")
dir.create(datadir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

baseline <- c(R0 = 2.5, K = 1e4, r = 0.125)
dt <- 0.1
t_max <- 200
nsim <- 1000L
n_patch <- 1L
alpha <- 0
gamma <- 1
I_init <- 10
tolerance <- 1e-10

required_files <- c(input_file, run_script, summary_script)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Required files not found: ", paste(missing_files, collapse = ", "))
}

grid <- readRDS(input_file)
required_columns <- c(
  "R0", "K", "r", "mean_ext_time.I1", "ext_prob.I1"
)
if (!is.data.frame(grid)) stop("Input RDS is not a data frame")
missing_columns <- setdiff(required_columns, names(grid))
if (length(missing_columns)) {
  stop("Input grid is missing columns: ", paste(missing_columns, collapse = ", "))
}

## Confirm simulation metadata from the executable array script.
run_text <- gsub(
  "[[:space:]]+", "",
  paste(readLines(run_script, warn = FALSE), collapse = "")
)
metadata_checks <- c(
  nsim = grepl("nsim<-1000L", run_text, fixed = TRUE),
  dt = grepl("dt<-if(reedfrost==1)1else0.1", run_text, fixed = TRUE),
  t_max = grepl("nt=round(200/dt)", run_text, fixed = TRUE),
  n_patch = grepl("n_patch=1", run_text, fixed = TRUE),
  alpha = grepl("alpha=0", run_text, fixed = TRUE),
  I_init = grepl("I_init=c(10,0)", run_text, fixed = TRUE),
  gamma = grepl("gamma=c(1,1)", run_text, fixed = TRUE),
  logistic = grepl("logistic_growth=logistic_growth", run_text, fixed = TRUE)
)
if (!all(metadata_checks)) {
  stop("Simulation metadata checks failed for: ",
       paste(names(metadata_checks)[!metadata_checks], collapse = ", "))
}

## sumfun_discrete() records the row index of first I1 == 0, averaged only
## over extinct runs. The first odin row is time zero, so row index j is
## converted to disease generations as (j - 1) * dt.
summary_text <- gsub(
  "[[:space:]]+", "",
  paste(readLines(summary_script, warn = FALSE), collapse = "")
)
conversion_verified <-
  grepl("first_ext<-function", summary_text, fixed = TRUE) &&
  grepl("idx[1L]", summary_text, fixed = TRUE) &&
  grepl(
    "mean_ext_time.I1=mean(ext1,na.rm=TRUE)",
    summary_text, fixed = TRUE
  )
if (!conversion_verified) {
  stop("Could not verify mean_ext_time.I1 row-index semantics")
}

near <- function(x, value) {
  abs(x - value) <= tolerance * max(1, abs(value))
}

expected_R0 <- seq(1.1, 5, by = 0.1)
expected_K <- 10^seq(3, 6, by = 0.25)
expected_r <- c(0.05, 0.1, 0.125, 0.2, 0.4)

same_numeric_set <- function(observed, expected) {
  observed <- sort(unique(observed))
  expected <- sort(expected)
  length(observed) == length(expected) &&
    all(abs(observed - expected) <=
          tolerance * pmax(1, abs(expected)))
}

grid_checks <- c(
  rows_2600 = nrow(grid) == 2600L,
  R0_values = same_numeric_set(grid$R0, expected_R0),
  K_values = same_numeric_set(grid$K, expected_K),
  r_values = same_numeric_set(grid$r, expected_r),
  no_duplicate_cells = !anyDuplicated(grid[c("R0", "K", "r")])
)
if (!all(grid_checks)) {
  stop("Grid checks failed for: ",
       paste(names(grid_checks)[!grid_checks], collapse = ", "))
}

clean_grid <- grid |>
  transmute(
    R0, K, log10K = log10(K), r,
    gamma = gamma, dt = dt, t_max = t_max, nsim = nsim,
    n_patch = n_patch, alpha = alpha, I_init = I_init,
    mean_ext_time_index = mean_ext_time.I1,
    mean_T_ext = (mean_ext_time.I1 - 1) * dt,
    ext_prob = ext_prob.I1,
    n_ext = round(nsim * ext_prob.I1),
    finite_conditional_mean =
      is.finite(mean_ext_time.I1) & (mean_ext_time.I1 - 1) * dt >= 0
  )

baseline_rows <- clean_grid |>
  filter(
    near(R0, baseline[["R0"]]),
    near(K, baseline[["K"]]),
    near(r, baseline[["r"]])
  )
if (nrow(baseline_rows) != 1L) {
  stop("Baseline R0=2.5, K=10000, r=0.125 is not uniquely present")
}

sweep_R0 <- clean_grid |>
  filter(near(K, baseline[["K"]]), near(r, baseline[["r"]])) |>
  arrange(R0)
sweep_K <- clean_grid |>
  filter(near(R0, baseline[["R0"]]), near(r, baseline[["r"]])) |>
  arrange(K)
sweep_r <- clean_grid |>
  filter(near(R0, baseline[["R0"]]), near(K, baseline[["K"]])) |>
  arrange(r)

sweep_checks <- c(
  R0_nonempty = nrow(sweep_R0) > 0L,
  K_nonempty = nrow(sweep_K) > 0L,
  r_nonempty = nrow(sweep_r) > 0L,
  R0_count = nrow(sweep_R0) == length(expected_R0),
  K_count = nrow(sweep_K) == length(expected_K),
  r_count = nrow(sweep_r) == length(expected_r),
  R0_fixed_values =
    all(near(sweep_R0$K, baseline[["K"]])) &&
    all(near(sweep_R0$r, baseline[["r"]])),
  K_fixed_values =
    all(near(sweep_K$R0, baseline[["R0"]])) &&
    all(near(sweep_K$r, baseline[["r"]])),
  r_fixed_values =
    all(near(sweep_r$R0, baseline[["R0"]])) &&
    all(near(sweep_r$K, baseline[["K"]])),
  extinction_probability_range =
    all(clean_grid$ext_prob >= 0 & clean_grid$ext_prob <= 1),
  conversion_finite_rows = all(
    abs(
      clean_grid$mean_T_ext[clean_grid$finite_conditional_mean] -
        (clean_grid$mean_ext_time_index[
          clean_grid$finite_conditional_mean
        ] - 1) * dt
    ) < 1e-12
  )
)
validation <- c(
  grid_checks,
  metadata_matches_source = all(metadata_checks),
  extinction_time_conversion_verified = conversion_verified,
  baseline_unique = nrow(baseline_rows) == 1L,
  sweep_checks
)
if (!all(validation)) {
  stop("Validation failed for: ",
       paste(names(validation)[!validation], collapse = ", "))
}

write.csv(
  clean_grid,
  file.path(datadir, "euler_onepatch_onestrain_extinct_demoggrid_clean.csv"),
  row.names = FALSE
)
write.csv(
  sweep_R0, file.path(datadir, "vary_R0_baseline.csv"), row.names = FALSE
)
write.csv(
  sweep_K, file.path(datadir, "vary_K_baseline.csv"), row.names = FALSE
)
write.csv(
  sweep_r, file.path(datadir, "vary_r_baseline.csv"), row.names = FALSE
)

make_ofat_figure <- function(data, x, x_label, varied_parameter,
                             fixed_text, baseline_x, x_scale = NULL) {
  p_time <- ggplot(data, aes(x = .data[[x]], y = mean_T_ext)) +
    geom_line(linewidth = 0.7, colour = "#0072B2", na.rm = TRUE) +
    geom_point(size = 2.1, colour = "#0072B2", na.rm = TRUE) +
    geom_vline(
      xintercept = baseline_x, colour = "grey35",
      linewidth = 0.5, linetype = "dashed"
    ) +
    labs(
      x = NULL,
      y = "Conditional mean extinction time\n(disease generations)",
      title = sprintf(
        "One-patch stochastic extinction: varying %s", varied_parameter
      ),
      subtitle = paste0(
        fixed_text,
        "; n_patch = 1; alpha = 0; gamma = 1; dt = 0.1; ",
        "t_max = 200; nsim = 1000; Poisson mean I(0) = 10"
      )
    ) +
    theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      plot.subtitle = element_text(size = 9),
      plot.margin = margin(8, 8, 4, 40, unit = "pt")
    )

  p_probability <- ggplot(data, aes(x = .data[[x]], y = ext_prob)) +
    geom_line(linewidth = 0.7, colour = "#D55E00") +
    geom_point(size = 2.1, colour = "#D55E00") +
    geom_vline(
      xintercept = baseline_x, colour = "grey35",
      linewidth = 0.5, linetype = "dashed"
    ) +
    scale_y_continuous(
      limits = c(0, 1), breaks = seq(0, 1, 0.2),
      labels = scales::percent_format(accuracy = 1)
    ) +
    labs(
      x = x_label,
      y = expression("Extinction probability by " * t[max] * " = 200"),
      caption = paste0(
        "Dashed vertical line marks the baseline ", varied_parameter,
        ". Mean extinction time is conditional on extinction by t_max."
      )
    ) +
    theme(
      plot.caption = element_text(hjust = 0, size = 8),
      plot.margin = margin(4, 8, 8, 40, unit = "pt")
    )

  if (!is.null(x_scale)) {
    p_time <- p_time + x_scale
    p_probability <- p_probability + x_scale
  }

  plot_grid(
    p_time, p_probability, ncol = 1, rel_heights = c(1.15, 1)
  )
}

figure_R0 <- make_ofat_figure(
  sweep_R0, "R0", expression(R[0]), "R0",
  "Fixed K = 10,000; r = 0.125",
  baseline[["R0"]]
)
figure_K <- make_ofat_figure(
  sweep_K, "K", "Carrying capacity K (log scale)", "K",
  "Fixed R0 = 2.5; r = 0.125",
  baseline[["K"]],
  scale_x_log10(
    breaks = 10^(3:6),
    labels = scales::label_number(big.mark = ",")
  )
)
figure_r <- make_ofat_figure(
  sweep_r, "r", "Host growth rate r", "r",
  "Fixed R0 = 2.5; K = 10,000",
  baseline[["r"]]
)

ggsave(
  file.path(figdir, "extinction_vary_R0.pdf"),
  figure_R0, width = 10, height = 8
)
ggsave(
  file.path(figdir, "extinction_vary_K.pdf"),
  figure_K, width = 10, height = 8
)
ggsave(
  file.path(figdir, "extinction_vary_r.pdf"),
  figure_r, width = 10, height = 8
)

cat("Validation checks passed: ", length(validation), "/",
    length(validation), "\n", sep = "")
cat("Baseline: R0 = 2.5, K = 10000, r = 0.125\n")
cat("Time conversion: mean_T_ext = (mean_ext_time.I1 - 1) * 0.1\n")
cat("Sweep points: R0 = ", nrow(sweep_R0), ", K = ", nrow(sweep_K),
    ", r = ", nrow(sweep_r), "\n", sep = "")
cat("No simulations were run. Outputs saved under: ", outdir, "\n", sep = "")
