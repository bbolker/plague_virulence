## Resumable stochastic validation runner.
## Examples:
##   Rscript fadeout/logistic_burnout/stochastic_validation/run_stochastic_validation.R --mode=main
##   Rscript fadeout/logistic_burnout/stochastic_validation/run_stochastic_validation.R --mode=r_sensitivity --min_sim=2000 --max_sim=10000

source(file.path(
  "fadeout", "logistic_burnout", "stochastic_validation",
  "stochastic_validation_functions.R"
))

module_dir <- file.path("fadeout", "logistic_burnout", "stochastic_validation")
output_dir <- ensure_dir(file.path(module_dir, "outputs"))

opt <- parse_cli_options(list(
  mode = "main",
  seed = 20260804,
  engine = "adaptive_tau",
  dt = 0.02,
  delta = 1e-6,
  initial_batch = 1000,
  min_sim = 5000,
  max_sim = 50000,
  target_ci_width = 0.01,
  tmax = 200,
  cores = 1,
  output_dir = file.path(module_dir, "outputs"),
  overwrite = FALSE
))
output_dir <- ensure_dir(as.character(opt$output_dir))

mode <- as.character(opt$mode)
engine <- as.character(opt$engine)
if (!engine %in% c("adaptive_tau", "tau_leap", "gillespie")) {
  stop("engine must be adaptive_tau, tau_leap, or gillespie")
}
grid <- make_parameter_grid(mode)
full_file <- file.path(output_dir, "full_stochastic_validation.csv")
boundary_file <- file.path(output_dir, "boundary_start_validation.csv")
counts_file <- file.path(output_dir, "outcome_counts.csv")
status_file <- file.path(output_dir, "validation_status_summary.csv")

done <- data.frame()
if (file.exists(full_file) && !as_logical(opt$overwrite)) {
  done <- read.csv(full_file)
}

full_rows <- if (nrow(done)) split(done, seq_len(nrow(done))) else list()
boundary_existing <- if (file.exists(boundary_file) && !as_logical(opt$overwrite)) {
  split(read.csv(boundary_file), seq_len(nrow(read.csv(boundary_file))))
} else {
  list()
}
boundary_rows <- boundary_existing

is_done <- function(p) {
  if (!nrow(done)) return(FALSE)
  any(done$R0 == p$R0 & done$r == p$r & done$K == p$K &
        done$I0 == p$I0 & done$delta == as_num(opt$delta) &
        done$engine == engine & done$dt == as_num(opt$dt))
}

for (i in seq_len(nrow(grid))) {
  p <- grid[i, ]
  if (is_done(p)) {
    message("Skipping completed cell ", i, "/", nrow(grid))
    next
  }
  seed <- cell_seed(as_int(opt$seed), p$R0, p$r, p$K, p$I0, i)
  message("Running cell ", i, "/", nrow(grid),
          ": R0=", p$R0, ", r=", p$r, ", K=", p$K)
  full <- adaptive_full_cell(
    p$R0, p$r, p$K, p$I0, delta = as_num(opt$delta),
    engine = engine, dt = as_num(opt$dt), seed = seed,
    initial_batch = as_int(opt$initial_batch),
    min_sim = as_int(opt$min_sim), max_sim = as_int(opt$max_sim),
    target_half_width = as_num(opt$target_ci_width) / 2,
    tmax = as_num(opt$tmax)
  )
  boundary <- simulate_tau_leap_boundary(
    p$R0, p$r, p$K, p$I0, n_sim = full$n_total,
    dt = as_num(opt$dt), seed = seed + 100000L, tmax = as_num(opt$tmax)
  )
  full_rows[[length(full_rows) + 1L]] <- full
  boundary_rows[[length(boundary_rows) + 1L]] <- boundary
  full_tab <- do.call(rbind, full_rows)
  boundary_tab <- do.call(rbind, boundary_rows)
  check_full_table(full_tab)
  write.csv(full_tab, full_file, row.names = FALSE)
  write.csv(boundary_tab, boundary_file, row.names = FALSE)
}

full_tab <- do.call(rbind, full_rows)
boundary_tab <- do.call(rbind, boundary_rows)
check_full_table(full_tab)
write.csv(full_tab, full_file, row.names = FALSE)
write.csv(boundary_tab, boundary_file, row.names = FALSE)

counts <- full_tab[c(
  "R0", "r", "K", "I0", "delta", "engine", "dt", "n_total",
  "n_fizzle", "n_late_extinction_before_boundary", "n_burnout",
  "n_persistence", "n_unresolved"
)]
write.csv(counts, counts_file, row.names = FALSE)

status <- data.frame(
  mode = mode,
  engine = engine,
  dt = as_num(opt$dt),
  delta = as_num(opt$delta),
  cells_completed = nrow(full_tab),
  total_realizations = sum(full_tab$n_total),
  unresolved_fraction_max = max(full_tab$unresolved_fraction, na.rm = TRUE),
  late_extinction_fraction_max = max(full_tab$late_extinction_fraction, na.rm = TRUE),
  mean_abs_error_unconditional = mean(full_tab$abs_error_unconditional, na.rm = TRUE),
  max_abs_error_unconditional = max(full_tab$abs_error_unconditional, na.rm = TRUE),
  stringsAsFactors = FALSE
)
metadata <- session_metadata()
status$git_sha <- metadata$git_sha[1]
write.csv(status, status_file, row.names = FALSE)

cat("Validation run completed\n")
print(status)
