## Quick local smoke test for the stochastic logistic-burnout validation.
## Run from repository root:
##   Rscript fadeout/logistic_burnout/stochastic_validation/run_smoke_test.R

source(file.path(
  "fadeout", "logistic_burnout", "stochastic_validation",
  "stochastic_validation_functions.R"
))

module_dir <- file.path("fadeout", "logistic_burnout", "stochastic_validation")
output_dir <- ensure_dir(file.path(module_dir, "outputs"))

opt <- parse_cli_options(list(
  seed = 20260804,
  n_sim = 150,
  dt = 0.02,
  delta = 1e-6,
  tmax = 200,
  cores = 1,
  output_dir = file.path(module_dir, "outputs")
))
output_dir <- ensure_dir(as.character(opt$output_dir))

grid <- make_parameter_grid("smoke")
rows <- vector("list", nrow(grid))
boundary_rows <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  p <- grid[i, ]
  seed <- cell_seed(as_int(opt$seed), p$R0, p$r, p$K, p$I0, i)
  message("Smoke full cell ", i, "/", nrow(grid),
          ": R0=", p$R0, ", r=", p$r, ", K=", p$K)
  rows[[i]] <- simulate_adaptive_tau_full(
    p$R0, p$r, p$K, p$I0, n_sim = as_int(opt$n_sim),
    delta = as_num(opt$delta), seed = seed, tmax = as_num(opt$tmax)
  )
  boundary_rows[[i]] <- simulate_tau_leap_boundary(
    p$R0, p$r, p$K, p$I0, n_sim = as_int(opt$n_sim),
    dt = as_num(opt$dt), seed = seed + 100000L, tmax = as_num(opt$tmax)
  )
}

full <- do.call(rbind, rows)
boundary <- do.call(rbind, boundary_rows)
check_full_table(full)

metadata <- session_metadata()
metadata$global_seed <- as_int(opt$seed)
metadata$mode <- "smoke"
metadata$dt <- as_num(opt$dt)
metadata$delta <- as_num(opt$delta)
metadata$n_sim <- as_int(opt$n_sim)

write.csv(full, file.path(output_dir, "smoke_test_results.csv"), row.names = FALSE)
write.csv(boundary, file.path(output_dir, "smoke_boundary_start_results.csv"),
          row.names = FALSE)
write.csv(metadata, file.path(output_dir, "smoke_metadata.csv"), row.names = FALSE)

status <- data.frame(
  check = c("probabilities_and_counts", "all_cells_completed"),
  passed = c(TRUE, nrow(full) == nrow(grid)),
  detail = c("check_full_table passed", paste(nrow(full), "cells"))
)
write.csv(status, file.path(output_dir, "smoke_status_summary.csv"),
          row.names = FALSE)

cat("Smoke test completed\n")
cat("Cells: ", nrow(full), "\n", sep = "")
cat("Total stochastic realizations: ", sum(full$n_total), "\n", sep = "")
print(full[c(
  "R0", "r", "K", "n_total", "n_fizzle",
  "n_late_extinction_before_boundary", "n_burnout",
  "n_persistence", "n_unresolved", "P_uncond_sim",
  "P_uncond_approx"
)])
