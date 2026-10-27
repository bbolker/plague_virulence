## Does the initial susceptible density affect single-strain persistence?
##
## Regime as in pairwise_invasion_100y.R (K = 5000, alpha = 5e-5/day,
## seasonal amplitude 0.2, 200 patches). Every patch starts at the same
## S = f * (S* + I*), capped at K, and I = I*, where (S*, I*) is the
## deterministic endemic equilibrium at the mean transmission rate, for f from
## 0.5 to 1.5. Strain 2 is inactive.
##
## Runs are simulated in 5-year chunks (every THIN-th day returned) and stop
## at global extinction, which is absorbing because immigration is
## proportional to total I. Persistence is global infection at 100 years.
##
## Usage:
##   Rscript fadeout/seasonal/s0_persistence_check.R
##
## Writes data/s0_persistence.rds (one row per run) and prints summaries;
## s0_persistence_check_plots.R draws the figure.

## Exec'd subprocesses inherit these; RhpcBLASctl caps this process
Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1, MKL_NUM_THREADS = 1)

library(dplyr)
library(tidyr)
library(purrr)
library(here)
library(odin)
library(future)
library(furrr)

## Each run gets its own L'Ecuyer-CMRG stream from furrr (seed = TRUE)
set.seed(20261015)

RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)


## ------------------------------------------------------------
## Parameters
## ------------------------------------------------------------

gamma <- 0.2              ## per day
r <- 0.125 * gamma        ## 0.125 per infectious period
n_patch <- 200
K <- 5000
alpha <- 5e-5             ## per day
seasonal_amp <- 0.2

years <- 100
chunk_days <- 5 * 365
THIN <- 10

R0_vals <- c(1.1, 1.2, 1.5, 2, 2.5, 3, 3.5, 5)
f_vals <- seq(0.5, 1.5, by = 0.1)
n_reps <- 20

## Jobs are spread over N_WORKERS forked workers (environment variable)
n_workers <- as.integer(Sys.getenv("N_WORKERS", 6))

grid <- expand_grid(R0 = R0_vals, f = f_vals, replicate = seq_len(n_reps))

outdir <- here::here("fadeout", "output", "s0_persistence_check")
dir.create(file.path(outdir, "data"), recursive = TRUE, showWarnings = FALSE)


## ------------------------------------------------------------
## Compile the shared odin model
## ------------------------------------------------------------

model_file <- here::here("fadeout", "seasonal", "seasonal_model_metapop.R")

## odin::odin() substitutes its argument, so the path must be a variable.
## The compiled model reaches the workers by forking, hence plan(multicore).
gen <- suppressMessages(odin::odin(model_file))

I_cols <- sprintf("I[%d,1]", seq_len(n_patch))


## ------------------------------------------------------------
## One realization
## ------------------------------------------------------------

run_one <- function(R0, f, replicate) {
  beta <- R0 * gamma
  S_star <- gamma * K / beta
  I_star <- r * (K - S_star) / beta
  S0 <- round(min(K, f * (S_star + I_star)))
  I0 <- round(I_star)
  mod <- gen$new(
    beta = c(beta, 0), gamma = c(gamma, gamma), dt = 1,
    I_ini = cbind(rep(I0, n_patch), 0), S_ini = rep(S0, n_patch),
    I2_ini = rep(0, n_patch), alpha = alpha,
    strain2_delay = .Machine$integer.max,
    r = rep(r, n_patch), K = rep(K, n_patch),
    season_period = 365, seasonal_amp = seasonal_amp, peak_day = 15,
    n_patch = n_patch
  )
  y <- NULL
  t0 <- 0L
  ext_day <- NA_real_
  occ <- numeric(0)
  ## Each chunk restarts from the previous chunk's final state
  while (t0 < years * 365) {
    steps <- seq(t0, t0 + chunk_days, by = THIN)
    raw <- mod$run(steps, y = y)
    glob <- rowSums(raw[, I_cols])
    occ <- c(occ, rowSums(raw[-1, I_cols] > 0) / n_patch)
    if (any(glob == 0)) {
      ext_day <- steps[match(0, glob)]
      break
    }
    y <- raw[nrow(raw), -1]
    t0 <- t0 + chunk_days
  }
  tibble(S0 = S0, I0 = I0, ext_day = ext_day,
         min_occ = min(occ), mean_occ = mean(occ))
}


## ------------------------------------------------------------
## Run, save, summarize
## ------------------------------------------------------------

plan(multicore(workers = n_workers))
res <- future_pmap(grid, run_one, .options = furrr_options(seed = TRUE))
plan(sequential)

out <- bind_cols(grid, list_rbind(res))
saveRDS(out, file.path(outdir, "data", "s0_persistence.rds"))

summ <- out |>
  group_by(R0, f, S0) |>
  summarise(p_persist = mean(is.na(ext_day)),
            med_ext_yr = median(ext_day / 365, na.rm = TRUE),
            .groups = "drop")

options(width = 200)
print(summ |> select(R0, f, p_persist) |>
        pivot_wider(names_from = f, values_from = p_persist))
print(summ |> select(R0, f, med_ext_yr) |>
        pivot_wider(names_from = f, values_from = med_ext_yr), digits = 2)

cat("Output: ", outdir, "\n", sep = "")
