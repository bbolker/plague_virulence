## Timestep and exact-engine checks for representative cells.
## Run from repository root:
##   Rscript fadeout/logistic_burnout/stochastic_validation/run_engine_checks.R

source(file.path(
  "fadeout", "logistic_burnout", "stochastic_validation",
  "stochastic_validation_functions.R"
))

module_dir <- file.path("fadeout", "logistic_burnout", "stochastic_validation")
output_dir <- ensure_dir(file.path(module_dir, "outputs"))
opt <- parse_cli_options(list(
  seed = 20260804,
  n_tau = 1000,
  n_gillespie = 400,
  delta = 1e-6,
  tmax = 200,
  cores = 1,
  output_dir = file.path(module_dir, "outputs")
))
output_dir <- ensure_dir(as.character(opt$output_dir))

grid <- make_parameter_grid("exact_subset")
dt_values <- c(0.05, 0.02, 0.01)
rows <- list()
idx <- 0L
for (i in seq_len(nrow(grid))) {
  p <- grid[i, ]
  for (dt in dt_values) {
    idx <- idx + 1L
    seed <- cell_seed(as_int(opt$seed), p$R0, p$r, p$K, p$I0, idx)
    message("Tau-leap engine check ", idx)
    rows[[length(rows) + 1L]] <- simulate_tau_leap_full(
      p$R0, p$r, p$K, p$I0, n_sim = as_int(opt$n_tau),
      dt = dt, delta = as_num(opt$delta), seed = seed,
      tmax = as_num(opt$tmax)
    )
  }
  idx <- idx + 1L
  seed <- cell_seed(as_int(opt$seed), p$R0, p$r, p$K, p$I0, idx)
  message("Adaptive tau engine check for R0=", p$R0, ", K=", p$K)
  rows[[length(rows) + 1L]] <- simulate_adaptive_tau_full(
    p$R0, p$r, p$K, p$I0, n_sim = as_int(opt$n_tau),
    delta = as_num(opt$delta), seed = seed, tmax = as_num(opt$tmax)
  )
  idx <- idx + 1L
  seed <- cell_seed(as_int(opt$seed), p$R0, p$r, p$K, p$I0, idx)
  message("Gillespie engine check for R0=", p$R0, ", K=", p$K)
  rows[[length(rows) + 1L]] <- simulate_gillespie_full(
    p$R0, p$r, p$K, p$I0, n_sim = as_int(opt$n_gillespie),
    delta = as_num(opt$delta), seed = seed, tmax = as_num(opt$tmax)
  )
}
tab <- do.call(rbind, rows)
check_full_table(tab)
write.csv(tab, file.path(output_dir, "simulation_engine_comparison.csv"),
          row.names = FALSE)
write.csv(tab[tab$engine == "tau_leap", ],
          file.path(output_dir, "dt_sensitivity.csv"), row.names = FALSE)
cat("Engine checks completed\n")
print(tab[c("R0", "r", "K", "engine", "dt", "n_total",
            "P_uncond_sim", "P_uncond_sim_ci_low",
            "P_uncond_sim_ci_high", "P_uncond_approx")])
