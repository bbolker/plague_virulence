## Burnout and persistence over a dense R0 x K x alpha grid, 100-year runs.
##
## Extends occupancy_factorial_100y.R to 8 K x 5 R0 x 7 alpha values, for
## seasonal (amplitude 0.40) and constant transmission, with replicate runs.
## Demography and initial conditions are as in that script: r = 0.125 per
## infectious period, gamma = 0.2/day, dt = 1 day, and each patch starts at a
## random point between a fresh outbreak and the endemic equilibrium.
##
## Patch-level diagnostics are computed inside each worker from the full daily
## state, because thinning hides short infection episodes. Only the saved
## occupancy time series is thinned (every THIN-th day).
##
## Episodes are maximal runs of days with I > 0 in a patch. Following
## extract_infection_episodes() in seasonal_fadeout_functions.R, an episode
## that ends within 1.5 T0 of its start (T0 = intrinsic period of the
## non-seasonal model) is an early extinction, and one that ends later is a
## fade-out. Early extinctions are split by whether the patch had a major
## outbreak (peak I >= major_frac * K):
##
##   failed                 early extinction, no major outbreak
##   burnout                early extinction after a major outbreak
##   fadeout                extinction after 1.5 T0
##   established_censored   still infected at the end, older than 1.5 T0
##   short_censored         still infected at the end, younger than 1.5 T0
##
## Usage:
##   Rscript fadeout/seasonal/occupancy_factorial_100y_batch.R [--pilot]
## --pilot runs only the corner combinations, one replicate each, and writes
## to a separate output directory.
##
## This script writes data only; occupancy_factorial_100y_batch_plots.R draws
## the heatmaps from data/combo_summary.csv and data/settings.csv.

library(plagueMetapop)
library(dplyr)
library(tidyr)
library(purrr)
library(here)
library(odin)
library(future)
library(furrr)

source(here::here("fadeout", "seasonal", "seasonal_fadeout_functions.R"))

pilot <- "--pilot" %in% commandArgs(trailingOnly = TRUE)

## Each run gets its own L'Ecuyer-CMRG stream from furrr (seed = TRUE)
set.seed(20261006)

RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)


## ------------------------------------------------------------
## Parameters
## ------------------------------------------------------------

dt <- 1                   ## days
gamma <- 0.2              ## per day
r <- 0.125 * gamma        ## 0.125 per infectious period
n_patch <- 200
season_period <- 365
peak_day <- 15

## Infected count of a freshly invaded patch in the ratchet initial conditions
I_outbreak <- 10

n_years <- 100
t_max <- n_years * 365
nt <- round(t_max / dt)

## Keep every THIN-th day of the saved occupancy series
THIN <- 10

## Episodes ending within threshold_multiplier * T0 are early extinctions
threshold_multiplier <- 1.5

## A major outbreak has peak I of at least major_frac * K
major_frac <- 0.01

## Late-run occupancy is summarized over the last summary_years years
summary_years <- 50

## Episode-duration histogram: bins of hist_width * T0, last bin is overflow
hist_width <- 0.1
hist_max <- 10

amps <- c(seasonal = 0.40, constant = 0.00)

n_reps <- 5

## Jobs are spread over at most this many forked workers, leaving at least
## one core free (detectCores() can return NA)
n_workers <- min(27, max(1, parallel::detectCores() - 1, na.rm = TRUE))

K_vals <- c(1000, 2000, 3000, 5000, 7500, 10000, 15000, 30000)
R0_vals <- seq(1.5, 3.5, by = 0.5)
alpha_vals <- c(1e-5, 2e-5, 3e-5, 5e-5, 8e-5, 1e-4, 2e-4)  ## per day

if (pilot) {
  K_vals <- range(K_vals)
  R0_vals <- range(R0_vals)
  alpha_vals <- range(alpha_vals)
  n_reps <- 1
}

param_grid <- expand_grid(
  R0 = R0_vals,
  K = K_vals,
  alpha = alpha_vals,
  scenario = names(amps),
  replicate = seq_len(n_reps)
) |>
  mutate(
    seasonal_amp = amps[scenario],
    run_id = row_number(),
    T0 = map_dbl(R0, \(x) calculate_intrinsic_period(x, gamma, r)$T0)
  )

cat("runs:", nrow(param_grid), "\n")


## ------------------------------------------------------------
## Output directory
## ------------------------------------------------------------

outdir <- here::here(
  "fadeout", "output",
  if (pilot) "occupancy_factorial_100y_batch_pilot"
  else "occupancy_factorial_100y_batch"
)

dir.create(file.path(outdir, "data"), recursive = TRUE, showWarnings = FALSE)


## ------------------------------------------------------------
## Compile the shared odin model
## ------------------------------------------------------------

model_file <- here::here("fadeout", "seasonal", "seasonal_model_metapop.R")

if (!file.exists(model_file)) {
  stop("Odin model file does not exist: ", model_file)
}

## odin::odin() substitutes its argument, so the path must be a variable.
## The compiled model is loaded into this process and reaches the workers by
## forking, which is why plan(multicore) is used rather than multisession.
gen <- suppressMessages(odin::odin(model_file))


## ------------------------------------------------------------
## Episode extraction and classification
## ------------------------------------------------------------

## I_mat: (nt + 1) x n_patch matrix of infected counts, row i is day i - 1
extract_episodes <- function(I_mat, T0, K) {
  n <- nrow(I_mat)
  threshold <- threshold_multiplier * T0
  eps <- lapply(seq_len(ncol(I_mat)), function(p) {
    Ip <- I_mat[, p]
    d <- diff(c(0L, as.integer(Ip > 0), 0L))
    starts <- which(d == 1L)
    ## first unoccupied row, or n + 1 if still infected at the end
    ends <- which(d == -1L)
    if (!length(starts)) return(NULL)
    stopifnot(sum(ends - starts) == sum(Ip > 0))
    data.frame(
      patch = p,
      start_day = starts - 1L,
      end_day = ifelse(ends > n, NA_integer_, ends - 1L),
      peak_I = vapply(seq_along(starts),
                      \(k) max(Ip[starts[k]:(ends[k] - 1L)]), numeric(1))
    )
  })
  bind_rows(eps) |>
    mutate(
      initial = start_day == 0L,
      censored = is.na(end_day),
      duration = ifelse(censored, nt - start_day, end_day - start_day),
      major = peak_I >= major_frac * K,
      class = case_when(
        !censored & duration <= threshold & !major ~ "failed",
        !censored & duration <= threshold ~ "burnout",
        !censored ~ "fadeout",
        duration > threshold ~ "established_censored",
        TRUE ~ "short_censored"
      )
    )
}

summarise_episodes <- function(eps, T0) {
  threshold <- threshold_multiplier * T0
  col <- eps |> filter(!initial)
  ini <- eps |> filter(initial)
  survived <- eps$duration > threshold
  tibble(
    n_episodes = nrow(eps),
    n_colonizations = nrow(col),
    n_failed = sum(col$class == "failed"),
    n_burnout = sum(col$class == "burnout"),
    n_survived = sum(col$duration > threshold),
    n_short_censored = sum(col$class == "short_censored"),
    ## fade-out hazard after establishment, from all episodes
    n_fadeout = sum(eps$class == "fadeout"),
    established_years = sum(eps$duration[survived] - threshold) / 365,
    n_local_extinctions = sum(!eps$censored),
    initial_failed = sum(ini$class == "failed"),
    initial_burnout = sum(ini$class == "burnout"),
    initial_survived = sum(ini$duration > threshold)
  ) |>
    mutate(
      p_failed = n_failed / (n_failed + n_burnout + n_survived),
      p_burnout = n_burnout / (n_burnout + n_survived),
      fadeout_hazard = n_fadeout / established_years
    )
}

duration_histogram <- function(eps, T0) {
  eps |>
    mutate(
      bin = pmin(floor(duration / T0 / hist_width), hist_max / hist_width)
    ) |>
    count(initial, class, bin, name = "n") |>
    mutate(bin = as.integer(bin), n = as.integer(n))
}


## ------------------------------------------------------------
## One realization
## ------------------------------------------------------------

run_one <- function(run_id, R0, K, alpha, seasonal_amp, T0, replicate, ...) {

  t0 <- Sys.time()
  beta0 <- R0 * gamma

  ## Deterministic equilibrium at the mean transmission rate
  S_star <- gamma * K / beta0
  I_star <- r * (K - S_star) / beta0

  ## Ratchet initial conditions: u = 1 is a fresh outbreak in a fully
  ## susceptible patch, u = 0 is the endemic equilibrium
  u <- runif(n_patch, min = 0, max = 1)
  S_ini <- round(K * u + S_star * (1 - u))
  I_ini <- cbind(round(I_outbreak * u + I_star * (1 - u)), rep(0, n_patch))

  mod <- gen$new(
    beta = c(beta0, 0),
    gamma = c(gamma, gamma),
    dt = dt,
    I_ini = I_ini,
    S_ini = S_ini,
    I2_ini = rep(0, n_patch),
    alpha = alpha,
    strain2_delay = .Machine$integer.max,
    r = rep(r, n_patch),
    K = rep(K, n_patch),
    season_period = season_period,
    seasonal_amp = seasonal_amp,
    peak_day = peak_day,
    n_patch = n_patch
  )

  raw <- mod$run(seq(0L, nt))
  I_mat <- raw[, sprintf("I[%d,1]", seq_len(n_patch))]
  S_mat <- raw[, sprintf("S[%d]", seq_len(n_patch))]
  rm(raw)

  occupied <- rowSums(I_mat > 0)
  global_I <- rowSums(I_mat)
  mean_N <- rowMeans(S_mat + I_mat)
  day <- seq(0L, nt)

  ## Immigration is proportional to total I, so global extinction is absorbing
  ext_row <- match(0, global_I)
  global_extinction_day <- if (is.na(ext_row)) NA_integer_ else day[ext_row]

  ## Patch-days before global extinction, for the Levins turnover rates
  at_risk_rows <- if (is.na(ext_row)) seq_along(day) else seq_len(ext_row - 1)
  occupied_years <- sum(occupied[at_risk_rows]) / 365
  empty_years <- (length(at_risk_rows) * n_patch) / 365 - occupied_years

  late <- day >= (n_years - summary_years) * 365
  late_occ <- occupied[late] / n_patch

  eps <- extract_episodes(I_mat, T0, K)
  ep_summary <- summarise_episodes(eps, T0)

  summary <- tibble(
    run_id = run_id,
    global_extinction_day = global_extinction_day,
    occ_mean = mean(late_occ),
    occ_median = median(late_occ),
    occ_lwr = unname(quantile(late_occ, 0.025)),
    occ_upr = unname(quantile(late_occ, 0.975)),
    final_occupied = occupied[length(occupied)],
    mean_N_late = mean(mean_N[late]),
    occupied_years = occupied_years,
    empty_years = empty_years
  ) |>
    bind_cols(ep_summary) |>
    mutate(
      ext_rate = n_local_extinctions / occupied_years,
      col_rate = n_colonizations / empty_years,
      elapsed_sec = as.numeric(difftime(Sys.time(), t0, units = "secs"))
    )

  keep <- seq(1, length(day), by = THIN)

  list(
    summary = summary,
    hist = duration_histogram(eps, T0) |> mutate(run_id = run_id, .before = 1),
    thinned = tibble(
      run_id = run_id,
      day = day[keep],
      occupied_patches = as.integer(occupied[keep]),
      global_I = as.integer(global_I[keep])
    ),
    episodes = if (replicate == 1) {
      eps |>
        mutate(run_id = run_id, .before = 1) |>
        mutate(across(c(start_day, end_day, duration, peak_I), as.integer))
    }
  )
}


## ------------------------------------------------------------
## Run all combinations in parallel
## ------------------------------------------------------------

plan(multicore(workers = min(n_workers, nrow(param_grid))))

t_start <- Sys.time()

res <- future_pmap(
  param_grid,
  run_one,
  .options = furrr_options(seed = TRUE)
)

plan(sequential)

total_elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
cat("total elapsed: ", round(total_elapsed, 1), " s\n", sep = "")

run_summary <- param_grid |>
  left_join(map(res, "summary") |> list_rbind(), by = "run_id")

duration_hist <- map(res, "hist") |> list_rbind()
occupancy_thinned <- map(res, "thinned") |> list_rbind()
episodes_rep1 <- map(res, "episodes") |> list_rbind()
rm(res)

cat("per-run elapsed (s):\n")
print(summary(run_summary$elapsed_sec))
cat("episodes per run:\n")
print(summary(run_summary$n_episodes))


## ------------------------------------------------------------
## Summaries across replicates
## ------------------------------------------------------------

## Probabilities and rates are pooled over replicates from summed counts and
## time at risk, not averaged per run
combo_summary <- run_summary |>
  group_by(R0, K, alpha, scenario, seasonal_amp, T0) |>
  summarise(
    n_runs = n(),
    p_persist = mean(is.na(global_extinction_day)),
    occ_median = median(occ_median),
    occ_mean = mean(occ_mean),
    p_failed = sum(n_failed) / sum(n_failed + n_burnout + n_survived),
    p_burnout = sum(n_burnout) / sum(n_burnout + n_survived),
    fadeout_hazard = sum(n_fadeout) / sum(established_years),
    ext_rate = sum(n_local_extinctions) / sum(occupied_years),
    col_rate = sum(n_colonizations) / sum(empty_years),
    p_initial_burnout = sum(initial_burnout) /
      sum(initial_burnout + initial_survived),
    .groups = "drop"
  )


## ------------------------------------------------------------
## Outputs
## ------------------------------------------------------------

write.csv(run_summary, file.path(outdir, "data", "run_summary.csv"),
          row.names = FALSE)
write.csv(combo_summary, file.path(outdir, "data", "combo_summary.csv"),
          row.names = FALSE)
saveRDS(duration_hist, file.path(outdir, "data", "duration_hist.rds"))
saveRDS(occupancy_thinned, file.path(outdir, "data", "occupancy_thinned.rds"))
saveRDS(episodes_rep1, file.path(outdir, "data", "episodes_rep1.rds"))

## Settings the plotting script needs for labels
settings <- mget(c(
  "n_patch", "n_years", "summary_years", "threshold_multiplier",
  "major_frac", "n_reps"
))
write.csv(
  tibble(name = names(settings), value = unlist(settings)),
  file.path(outdir, "data", "settings.csv"),
  row.names = FALSE
)

cat("Output: ", outdir, "\n", sep = "")
