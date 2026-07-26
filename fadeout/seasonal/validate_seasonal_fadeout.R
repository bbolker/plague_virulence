## Focused synthetic validation plus an optional reference-scenario smoke test.
## Run from the repository root:
##   Rscript fadeout/seasonal/validate_seasonal_fadeout.R
##   Rscript fadeout/seasonal/validate_seasonal_fadeout.R --smoke

library(here)
source(here::here(
  "fadeout", "seasonal", "seasonal_fadeout_functions.R"
))

smoke <- "--smoke" %in% commandArgs(trailingOnly = TRUE)

period <- calculate_intrinsic_period(R0 = 2.5, gamma = 0.2, r = 0.02)
stopifnot(abs(period$omega0 - period$eigen_omega0) < 1e-10)
threshold_multiplier <- 1.5
threshold <- threshold_multiplier * period$T0
dt <- 10
times <- seq(0, 500, by = dt)

## Four patches encode: early extinction, later fade-out, established censoring,
## and two distinct episodes separated by an observed zero interval.
synthetic <- rbind(
  data.frame(
    time = times, patch = 1L,
    I = as.integer(times >= 0 & times < floor(threshold / dt) * dt)
  ),
  data.frame(
    time = times, patch = 2L,
    I = as.integer(times >= 0 & times < ceiling((threshold + 30) / dt) * dt)
  ),
  data.frame(
    time = times, patch = 3L,
    I = as.integer(times >= 250)
  ),
  data.frame(
    time = times, patch = 4L,
    I = as.integer((times >= 0 & times < 20) | (times >= 100 & times < 120))
  )
)
episodes <- extract_infection_episodes(
  synthetic, "synthetic", 1L, period$T0, threshold_multiplier
)

patch1 <- episodes[episodes$patch == 1, ]
patch2 <- episodes[episodes$patch == 2, ]
patch3 <- episodes[episodes$patch == 3, ]
patch4 <- episodes[episodes$patch == 4, ]
stopifnot(
  nrow(patch1) == 1L, patch1$early_extinction, !patch1$fadeout,
  nrow(patch2) == 1L, patch2$survived_beyond_threshold, patch2$fadeout,
  nrow(patch3) == 1L, patch3$survived_beyond_threshold,
  patch3$right_censored, !patch3$fadeout,
  nrow(patch4) == 2L
)

occupancy_known <- data.frame(
  time = 0:8,
  occupied_patches = c(0, 1, 2, 1, 2, 2, 0, 0, 0),
  occupancy = c(0, 0.5, 1, 0.5, 1, 1, 0, 0, 0),
  total_I = c(0, 1, 2, 1, 2, 2, 0, 0, 0)
)
annual <- summarise_annual_occupancy(
  occupancy_known, "known", 1L, season_period = 4, dt = 1
)
stopifnot(
  abs(annual$mean_occupancy[1] - 0.5) < 1e-12,
  abs(annual$occupancy_amplitude[1] - 1) < 1e-12,
  abs(annual$mean_occupancy[2] - 0.4) < 1e-12,
  abs(annual$occupancy_amplitude[2] - 1) < 1e-12
)
cat("Synthetic validation passed.\n")

if (smoke) {
  smoke_grid <- data.frame(
    combo_id = "smoke", R0 = 2.5, gamma = 0.2, r = 0.02,
    K = 3000, alpha = 1e-4, seasonal_amp = 0.4,
    season_period = 365, peak_day = 15, n_patch = 200L, dt = 1,
    t_max = 10 * 365, n_reps = 1L, base_seed = 101L,
    threshold_multiplier = 1.5
  )
  output_dir <- here::here(
    "fadeout", "output", "seasonal_fadeout", "smoke"
  )
  run_seasonal_fadeout_grid(
    smoke_grid, output_dir,
    here::here("fadeout", "seasonal", "seasonal_model_metapop.R")
  )
  expected <- c(
    "data/infection_episodes.csv",
    "data/annual_occupancy.csv",
    "data/annual_fadeout_counts.csv",
    "data/replicate_summary.csv",
    "data/parameters_and_seeds.csv",
    "data/global_extinction.csv",
    "figures/annual_mean_occupancy.pdf",
    "figures/annual_occupancy_amplitude.pdf",
    "figures/episode_survival.pdf",
    "figures/established_episode_survival.pdf",
    "figures/occupancy_raster_smoke_rep001.pdf"
  )
  missing <- expected[!file.exists(file.path(output_dir, expected))]
  if (length(missing)) {
    stop("Smoke test did not create: ", paste(missing, collapse = ", "))
  }
  cat("Smoke-test simulation and output checks passed.\n")
  cat("Smoke output: ", output_dir, "\n", sep = "")
}
