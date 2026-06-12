## Run singlepatchintro simulations for the six parameter sets in
## euler_twostrain_example_pars.csv, using the same simulation settings as the
## SLURM job (sharcnet/euler_twostrain_examples_run_array.R).
## Per-step summaries computed via Arrow and stored as a list.

library(plagueMetapop)
library(future)
library(dplyr)
library(arrow)

## --- simulation constants (match euler_twostrain_run_array.R singlepatchintro)
n_patch       <- 200
n_sim         <- 500L
dt            <- 0.1
strain2_delay <- round(100 / dt)
nt            <- round(500 / dt)
stop_cond     <- stop_both_extinct
I_init        <- cbind(rep(10, n_patch),
                       c(10, rep(0, n_patch - 1L)))

pars <- read.csv(here::here("metapop/odin/euler_twostrain_example_pars.csv"))

summarise_runs <- function(runs_chunk) {
  lapply(runs_chunk, function(x) {
    x |>
      arrow_table() |>
      filter(step > strain2_delay * dt,
             abs(step %% 1.0) < 1e-6) |>
      summarise(patches = sum(value > 0),
                mean    = mean(if_else(value == 0, NA_real_, value), na.rm = TRUE),
                sd      = sd(if_else(value == 0, NA_real_, value),   na.rm = TRUE),
                .by = c(state, step)) |>
      collect()
  })
}

plan(multicore(workers = 4))

results <- vector("list", nrow(pars))

for (i in seq_len(nrow(pars))) {
  cat(i, "\n")
  row <- pars[i, ]
  print(row)
  chunk_size <- 25L
  n_chunks   <- n_sim %/% chunk_size
  chunks     <- vector("list", n_chunks)

  for (k in seq_len(n_chunks)) {
    runs_k    <- discrete_run(beta_vec      = c(row$R01, row$R02),
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
                              nsim          = chunk_size,
                              platform      = "odin")
    chunks[[k]] <- bind_rows(summarise_runs(runs_k), .id = "run") |>
      mutate(run = as.integer(run) + (k - 1L) * chunk_size)
    rm(runs_k); gc()
  }

  results[[i]] <- bind_rows(chunks)
}

plan(sequential)

attr(results, "metadata") <- pars

saveRDS(results,
        here::here("metapop/odin/outputs/euler_twostrain_singlepatchintro_examples.rds"))
