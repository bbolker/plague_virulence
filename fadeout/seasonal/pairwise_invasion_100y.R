## Pairwise invasion runs in a regime with intermediate long-term occupancy:
## K = 5000, seasonal amplitude 0.2, alpha = 5e-5 per day (see the
## occupancy_factorial_100y_batch heatmaps).
##
## For each (resident R0, invader R0) pair, the resident (strain 1) runs alone
## for 100 years from the ratchet initial conditions of
## occupancy_factorial_100y_batch.R. Then 10 infecteds of the invader
## (strain 2) are put into one patch, chosen at random among the patches the
## resident occupies (any patch if the resident is extinct), and both strains
## are followed for another 50 years. The simulation code and fixed settings
## are in pairwise_invasion_functions.R.
##
## Demography is as in occupancy_factorial_100y_batch.R: r = 0.125 per
## infectious period, gamma = 0.2/day, dt = 1 day, 200 patches.
##
## Outcomes at the end of the invasion phase:
##
##   resident_absent    resident globally extinct before the invasion
##   invader_lost       only the resident is left
##   resident_replaced  only the invader is left
##   coexist            both strains are present
##   both_extinct       neither strain is left
##
## Immigration is proportional to total I of each strain, so global extinction
## of a strain is absorbing.
##
## Usage:
##   Rscript fadeout/seasonal/pairwise_invasion_100y.R [--pilot] [--add-missing]
##     [--susc-eq]
## --pilot runs only the corner pairs, one replicate each, and writes to a
## separate output directory.
## --susc-eq starts S in each patch at an independent Uniform(0.7, 1.3)
## multiple of the endemic equilibrium population size N* = S* + I*, capped
## at K, instead of the ratchet S (initial I is unchanged), and adds
## "_susc_eq" to the output directory name.
## --add-missing reads the existing output, runs only the (R0_res, R0_inv,
## replicate) combinations it lacks (e.g. after R0 values are added), and
## writes the combined output. Without it, the whole grid is run. Combined
## output is not bit-identical to a run from scratch, because furrr's seed
## streams follow position in the grid of runs actually done.
##
## The factorial over K, alpha and amplitude runs on Nibi:
## odin/sharcnet/pairwise_invasion_factorial_run_array.R.
##
## This script writes data only; pairwise_invasion_100y_plots.R draws the
## figures.

## Exec'd subprocesses inherit these; RhpcBLASctl caps this process
Sys.setenv(OPENBLAS_NUM_THREADS = 1, OMP_NUM_THREADS = 1, MKL_NUM_THREADS = 1)

library(dplyr)
library(tidyr)
library(purrr)
library(here)
library(odin)
library(future)
library(furrr)

source(here::here("fadeout", "seasonal", "pairwise_invasion_functions.R"))

pilot <- "--pilot" %in% commandArgs(trailingOnly = TRUE)
add_missing <- "--add-missing" %in% commandArgs(trailingOnly = TRUE)
susc_eq <- "--susc-eq" %in% commandArgs(trailingOnly = TRUE)

## Each run gets its own L'Ecuyer-CMRG stream from furrr (seed = TRUE)
set.seed(20261015)

RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)


## ------------------------------------------------------------
## Parameters (fixed settings are in pairwise_invasion_functions.R)
## ------------------------------------------------------------

K <- 5000
alpha <- 5e-5             ## per day
seasonal_amp <- 0.2

R0_vals <- c(1.1, 1.2, seq(1.5, 3.5, by = 0.5), 5)
n_reps <- 100

## Jobs are spread over N_WORKERS forked workers (environment variable)
n_workers <- as.integer(Sys.getenv("N_WORKERS", 15))

if (pilot) {
  R0_vals <- range(R0_vals)
  n_reps <- 1
}

param_grid <- expand_grid(
  R0_res = R0_vals,
  R0_inv = R0_vals,
  replicate = seq_len(n_reps)
) |>
  mutate(run_id = row_number())

grid_keys <- c("R0_res", "R0_inv", "replicate")


## ------------------------------------------------------------
## Output directory, and existing runs for --add-missing
## ------------------------------------------------------------

outdir <- here::here(
  "fadeout", "output",
  paste0(
    if (pilot) "pairwise_invasion_100y_pilot" else "pairwise_invasion_100y",
    if (susc_eq) "_susc_eq"
  )
)

dir.create(file.path(outdir, "data"), recursive = TRUE, showWarnings = FALSE)

run_grid <- param_grid

if (add_missing) {
  old_summary <- read.csv(file.path(outdir, "data", "run_summary.csv"))
  old_thinned <- readRDS(file.path(outdir, "data", "occupancy_thinned.rds"))

  ## Map old run ids onto the current grid; drop runs no longer in it
  id_map <- old_summary |>
    select(all_of(grid_keys), old_id = run_id) |>
    inner_join(param_grid |> select(all_of(grid_keys), run_id),
               by = grid_keys)
  old_summary <- old_summary |>
    select(-run_id) |>
    inner_join(param_grid, by = grid_keys) |>
    relocate(names(param_grid))
  old_thinned <- old_thinned |>
    inner_join(id_map |> select(old_id, new_id = run_id),
               by = c(run_id = "old_id")) |>
    mutate(run_id = new_id, .keep = "unused")

  run_grid <- anti_join(param_grid, id_map, by = grid_keys)
  cat("existing runs kept:", nrow(old_summary), "\n")
}

cat("runs:", nrow(run_grid), "\n")


## ------------------------------------------------------------
## Run all pairs in parallel
## ------------------------------------------------------------

gen <- make_generator(
  here::here("fadeout", "seasonal", "seasonal_model_metapop.R")
)

plan(multicore(workers = max(1, min(n_workers, nrow(run_grid)))))

t_start <- Sys.time()

res <- future_pmap(
  run_grid,
  run_one,
  K = K, alpha = alpha, seasonal_amp = seasonal_amp, gen = gen,
  susc_eq = susc_eq,
  .options = furrr_options(seed = TRUE, scheduling = 4)
)

plan(sequential)

total_elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
cat("total elapsed: ", round(total_elapsed, 1), " s\n", sep = "")

run_summary <- run_grid |>
  inner_join(map(res, "summary") |> list_rbind(), by = "run_id")

occupancy_thinned <- map(res, "thinned") |> list_rbind()
rm(res)

cat("per-run elapsed (s), new runs:\n")
print(summary(run_summary$elapsed_sec))

if (add_missing) {
  run_summary <- bind_rows(old_summary, run_summary) |> arrange(run_id)
  occupancy_thinned <- bind_rows(old_thinned, occupancy_thinned) |>
    arrange(run_id, day)
  rm(old_summary, old_thinned)
}
stopifnot(nrow(run_summary) == nrow(param_grid))

cat("outcomes:\n")
print(table(run_summary$outcome))

combo_summary <- summarise_pairs(run_summary, R0_res, R0_inv)


## ------------------------------------------------------------
## Outputs
## ------------------------------------------------------------

write.csv(run_summary, file.path(outdir, "data", "run_summary.csv"),
          row.names = FALSE)
write.csv(combo_summary, file.path(outdir, "data", "combo_summary.csv"),
          row.names = FALSE)
saveRDS(occupancy_thinned, file.path(outdir, "data", "occupancy_thinned.rds"))

## Settings the plotting script needs for labels
settings <- mget(c(
  "n_patch", "K", "alpha", "seasonal_amp", "susc_eq", "resident_years",
  "invasion_years", "summary_years", "I_invader", "n_reps", "THIN"
))
write.csv(
  tibble(name = names(settings), value = unlist(settings)),
  file.path(outdir, "data", "settings.csv"),
  row.names = FALSE
)

cat("Output: ", outdir, "\n", sep = "")
