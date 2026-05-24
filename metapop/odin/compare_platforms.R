## Compare odin, macpan2, and pureR implementations of the two-strain n-patch model.
## Tests summary statistics (mean patch trajectories) rather than exact equality
## due to RNG differences between engines.
## For odin and macpan2, simulator creation and trajectory generation are timed separately.

suppressPackageStartupMessages({
  library(odin); library(dde); library(tidyr); library(dplyr)
  library(macpan2); library(ggplot2)
})
options(macpan2_verbose = FALSE)

source(here::here("metapop/odin", "discrete_odin.R"))
source(here::here("metapop/odin", "discrete_macpan2.R"))
source(here::here("metapop/odin", "discrete_pureR.R"))

## -- timings -------------------------------------------------------------------
cat("=== Timings (nt=1000, n_patch=100) ===\n")
set.seed(42)
t_odin_setup    <- system.time(sim_odin <- make_simulator_odin())
t_odin_traj     <- system.time(res_odin <- run_simulator_odin(sim_odin))
t_macpan2_setup <- system.time(sim_macpan2 <- make_simulator_macpan2(seed = 42))
t_macpan2_traj  <- system.time(res_macpan2 <- run_simulator_macpan2(sim_macpan2))
t_pureR         <- system.time(res_pureR   <- run_simulator_pureR(seed = 42))

cat(sprintf("odin    (setup):   %.2f s\n", t_odin_setup["elapsed"]))
cat(sprintf("odin    (traj):    %.2f s\n", t_odin_traj["elapsed"]))
cat(sprintf("macpan2 (setup):   %.2f s\n", t_macpan2_setup["elapsed"]))
cat(sprintf("macpan2 (traj):    %.2f s\n", t_macpan2_traj["elapsed"]))
cat(sprintf("pureR   (total):   %.2f s\n", t_pureR["elapsed"]))

## -- convert to common long format ---------------------------------------------
long_odin    <- conv_odin(res_odin)
long_macpan2 <- conv_macpan2(res_macpan2)
long_pureR   <- res_pureR  ## already in long format

## -- summaries: mean over patches at each time step ---------------------------
means_from <- function(df) df |>
  group_by(step, state) |>
  summarise(mean_val = mean(value, na.rm = TRUE), .groups = "drop")

odin_means    <- means_from(long_odin)
macpan2_means <- means_from(long_macpan2)
pureR_means   <- means_from(long_pureR)

fmt <- function(df) df |>
  group_by(state) |>
  summarise(grand_mean = round(mean(mean_val), 1),
            grand_sd   = round(sd(mean_val),   1))

cat("\n=== Grand mean +/- sd of patch-mean trajectory (all steps) ===\n")
cat("odin:\n");    print(fmt(odin_means))
cat("macpan2:\n"); print(fmt(macpan2_means))
cat("pureR:\n");   print(fmt(pureR_means))

cat("\n=== Patch mean at final time step ===\n")
cat("odin:\n");    print(odin_means    |> filter(step == max(step)))
cat("macpan2:\n"); print(macpan2_means |> filter(step == max(step)))
cat("pureR:\n");   print(pureR_means   |> filter(step == max(step)))
