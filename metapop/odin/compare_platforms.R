## Compare odin, macpan2, and pureR implementations of the two-strain n-patch model.
## Tests summary statistics (mean patch trajectories) rather than exact equality
## due to RNG differences between engines.
## For odin and macpan2, simulator creation and trajectory generation are timed separately.

suppressPackageStartupMessages({
  library(odin); library(dde); library(tidyr); library(dplyr)
  library(macpan2); library(ggplot2)
})
options(macpan2_verbose = FALSE)

platforms <- c("odin", "macpan2", "pureR")
for (p in platforms)  {
  source(here::here("metapop/odin", sprintf("discrete_%s.R", p)))
}

## FIXME: set up basic S3 class structure to avoid all the mget() nonsense


## -- timings -------------------------------------------------------------------
cat("=== Timings (nt=1000, n_patch=100) ===\n")
set.seed(42)
make_funs <- mget(sprintf("make_simulator_%s", platforms))
t_setup <- list()
sim_list <- list()
for (i in seq_along(make_funs)) {
  t_setup[[platforms[[i]]]] <- system.time(sim_list[[platforms[[i]]]] <- make_funs[[i]]())
}
t_traj <- list()
traj_list <- list()
for (p in platforms) {
  run_fun <- get(sprintf("run_simulator_%s", p))
  t_traj[[p]] <- system.time(traj_list[[p]] <- run_fun(sim_list[[p]]))
}

sapply(t_traj, \(x) x[["elapsed"]])
sapply(t_setup, \(x) x[["elapsed"]])

## -- convert to common long format ---------------------------------------------
conv_list <- list()
for (p in platforms) {
  conv_list[[p]] <- get(sprintf("conv_%s", p))(traj_list[[p]])
}

## mean across patches
mean_df <- conv_list |>
  bind_rows(.id = "platform") |>
  summarise(mean_val = mean(value, na.rm = TRUE),
            .by = c(platform, step, state))

## grand mean
mean_df |> summarise(across(mean_val, .fns = c(mean = mean, sd = sd), .names = "grand_{.fn}"),
                     .by = c(platform, state))

## last value
mean_df |> summarise(across(mean_val, ~tail(., 1)),
                     .by = c(platform, state))
