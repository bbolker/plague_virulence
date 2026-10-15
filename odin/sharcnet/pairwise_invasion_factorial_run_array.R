## HPC array version of fadeout/seasonal/pairwise_invasion_100y.R, over a
## K x alpha x seasonal amplitude factorial.
## Submit via submit_pairwise_invasion_factorial.sh (full, through
## submit_chunked.sh) or submit_pairwise_invasion_factorial_mini.sh (mini).
## Combine outputs with pairwise_invasion_factorial_combine.R
##
## Each task is one (K, alpha, seasonal_amp, R0_res, R0_inv) row of the grid:
## n_reps runs of a 100-year resident phase, a single-patch invasion, and 50
## years of two-strain dynamics (see fadeout/seasonal/pairwise_invasion_functions.R).
## The thinned occupancy series is kept only for the first n_save_series
## replicates.

library(dplyr)
library(purrr)
library(here)
library(future)
library(furrr)
library(optparse)

opt <- parse_args(OptionParser(option_list = list(
  make_option(c("-m", "--mini"), action = "store_true", default = FALSE,
              help = "run mini test grid (8 environments x 4 pairs, 5 reps)")
)))

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(task_id)) stop("SLURM_ARRAY_TASK_ID not set")

source(here::here("fadeout", "seasonal", "pairwise_invasion_functions.R"))

base_fn <- "pairwise_factorial"
if (opt$mini) base_fn <- paste0(base_fn, "_mini")

if (opt$mini) {
  dd <- expand.grid(R0_res       = c(1.5, 3),
                    R0_inv       = c(1.5, 3),
                    K            = c(5000, 10000),
                    alpha        = c(1e-5, 2e-4),
                    seasonal_amp = c(0.1, 0.4))
  n_reps <- 5L
} else {
  R0_vals <- c(1.2, 1.5, 2, 2.5, 3, 3.5, 5)
  dd <- expand.grid(R0_res       = R0_vals,
                    R0_inv       = R0_vals,
                    K            = c(5000, 7500, 10000),
                    alpha        = c(1e-5, 5e-5, 2e-4),
                    seasonal_amp = c(0.1, 0.2, 0.4))
  n_reps <- 100L
}
n_save_series <- 20L

if (task_id > nrow(dd)) stop("task ", task_id, " > grid size ", nrow(dd))
row <- dd[task_id, ]

out_file <- sprintf("outputs/%s_task_%06d.rds", base_fn, task_id)
if (file.exists(out_file)) {
  cat("output exists, skipping:", out_file, "\n")
  quit(save = "no")
}

## Each replicate gets its own L'Ecuyer-CMRG stream from furrr (seed = TRUE)
set.seed(20261016 + task_id)

gen <- make_generator(
  here::here("fadeout", "seasonal", "seasonal_model_metapop.R")
)

t_start <- Sys.time()
plan(multicore(workers = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 1L))))
res <- future_map(
  seq_len(n_reps),
  \(i) run_one(run_id = i, R0_res = row$R0_res, R0_inv = row$R0_inv,
               K = row$K, alpha = row$alpha,
               seasonal_amp = row$seasonal_amp, gen = gen),
  .options = furrr_options(seed = TRUE)
)
plan(sequential)
cat("task", task_id, "elapsed:",
    round(as.numeric(difftime(Sys.time(), t_start, units = "secs"))), "s\n")

## run_id within a task is the replicate number
summary <- map(res, "summary") |>
  list_rbind() |>
  rename(replicate = run_id)
thinned <- map(res[seq_len(min(n_save_series, n_reps))], "thinned") |>
  list_rbind() |>
  rename(replicate = run_id)

out <- list(
  row = bind_cols(task_id = task_id, row),
  summary = bind_cols(task_id = task_id, row, summary),
  thinned = bind_cols(task_id = task_id, thinned)
)

dir.create("outputs", showWarnings = FALSE)
saveRDS(out, out_file)
