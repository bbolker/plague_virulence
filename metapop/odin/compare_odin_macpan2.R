## Compare odin and macpan2 implementations of the two-strain n-patch model.
## Tests summary statistics (mean patch trajectories) rather than exact equality
## due to RNG differences between engines.

suppressPackageStartupMessages({
  library(odin); library(dde); library(tidyr); library(dplyr)
  library(macpan2); library(ggplot2)
})
options(macpan2_verbose = FALSE)

## -- odin setup ----------------------------------------------------------------
odin_file <- here::here("metapop/odin", "odin_twostrain_npatch.R")
gen_odin  <- odin::odin(odin_file)

run_odin <- function(beta_vec=c(1.5,2.5), K=1e4, r=0.125, n_patch=100,
                     n_strain=2, nt=1000, alpha=1e-3, I_init=10, seed=NULL) {
  K_vec   <- rep(K, length.out = n_patch)
  r_vec   <- rep(r, length.out = n_patch)
  I_init2 <- rep(I_init, length.out = 2)
  if (!is.null(seed)) set.seed(seed)
  I_ini_mat <- matrix(rpois(n_strain*n_patch, lambda=I_init2), byrow=TRUE, ncol=n_strain)
  S_ini_vec <- K_vec - rowSums(I_ini_mat)
  mod <- gen_odin$new(beta=beta_vec, r=r_vec, K=K_vec, I_ini=I_ini_mat,
                      S_ini=S_ini_vec, alpha=alpha, n_patch=n_patch)
  mod$run(1:nt)
}

conv_fun <- function(x) {
  as.data.frame.table(x[, !colnames(x) %in% "step"]) |>
    as_tibble() |>
    transmute(step  = as.numeric(Var1),
              state = case_when(
                grepl("^S", Var2) ~ gsub("\\[", ".", gsub("\\]", "", Var2)),
                grepl("^I", Var2) ~ gsub("\\[([0-9]+),([0-9]+)\\]", "\\2.\\1", Var2)),
              value = Freq) |>
    mutate(patch = as.numeric(sub(".*\\.", "", state)),
           state = sub("\\..*", "", state)) |>
    arrange(step, patch)
}

## -- macpan2 setup -------------------------------------------------------------
source(here::here("metapop/odin", "macpan2_twostrain_npatch.R"))

## -- timings -------------------------------------------------------------------
cat("=== Timings (nt=1000, n_patch=100) ===\n")
t_odin    <- system.time(res_odin    <- run_odin(seed = 42))
t_macpan2 <- system.time(res_macpan2 <- run_twostrain_macpan2(seed = 42))
cat(sprintf("odin:    %.2f s\n", t_odin["elapsed"]))
cat(sprintf("macpan2: %.2f s\n", t_macpan2["elapsed"]))

## -- summaries: mean over patches at each time step ---------------------------
odin_means <- conv_fun(res_odin) |>
  group_by(step, state) |>
  summarise(mean_val = mean(value, na.rm=TRUE), .groups="drop")

macpan2_means <- res_macpan2 |>
  group_by(step, state) |>
  summarise(mean_val = mean(value, na.rm=TRUE), .groups="drop")

fmt <- function(df) df |>
  group_by(state) |>
  summarise(grand_mean = round(mean(mean_val), 1),
            grand_sd   = round(sd(mean_val),   1))

cat("\n=== Grand mean +/- sd of patch-mean trajectory (all steps) ===\n")
cat("odin:\n");    print(fmt(odin_means))
cat("macpan2:\n"); print(fmt(macpan2_means))

cat("\n=== Patch mean at final time step ===\n")
cat("odin:\n");    print(odin_means    |> filter(step == max(step)))
cat("macpan2:\n"); print(macpan2_means |> filter(step == max(step)))
