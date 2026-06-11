## SLURM array version of euler_twostrain_singlepatchintro_examples.R
## One task per row of euler_twostrain_example_pars.csv (6 rows → array 1-6).
## Submit via submit_euler_twostrain_examples.sh.
## Combine outputs with euler_twostrain_examples_combine.R.

library(plagueMetapop)
library(future)
library(dplyr)

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(task_id)) stop("SLURM_ARRAY_TASK_ID not set")

pars <- read.csv("euler_twostrain_example_pars.csv")
row  <- pars[task_id, ]

## simulation constants (match euler_twostrain_run_array.R singlepatchintro)
n_patch       <- 200
n_sim         <- 500L
dt            <- 0.1
strain2_delay <- round(100 / dt)
nt            <- round(500 / dt)
stop_cond     <- stop_both_extinct
I_init        <- cbind(rep(10, n_patch),
                       c(10, rep(0, n_patch - 1L)))

plan(multicore(workers = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 1L))))
runs <- discrete_run(beta_vec      = c(row$R01, row$R02),
                     K             = row$K,
                     r             = 0.125,
                     n_patch       = n_patch,
                     nt            = nt,
                     alpha         = row$alpha,
                     I_init        = I_init,
                     gamma         = c(1, 1),
                     dt            = dt,
                     def_file      = "euler_odin_def.R",
                     strain2_delay = strain2_delay,
                     stop_cond     = stop_cond,
                     nsim          = n_sim,
                     platform      = "odin")
plan(sequential)

## thin each run before binding to keep peak memory low
result <- bind_rows(
  lapply(seq_along(runs), \(j) {
    thin_idx <- seq(1L, nrow(runs[[j]]), by = 10L)
    runs[[j]][thin_idx, ] |> mutate(run = as.character(j))
  })
) |> mutate(K = row$K, alpha = row$alpha, R01 = row$R01, R02 = row$R02,
            description = row$description)

dir.create("outputs", showWarnings = FALSE)
saveRDS(result, sprintf("outputs/euler_twostrain_examples_task_%06d.rds", task_id))
