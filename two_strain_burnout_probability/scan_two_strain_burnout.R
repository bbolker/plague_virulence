## Simple two-strain burnout parameter scan.
## Put this file in the same folder as two_strain_burnout.R, then run:
## source("scan_two_strain_burnout_v4.R")

if (!requireNamespace("burnout", quietly = TRUE)) {
  stop("Please install the burnout package: remotes::install_github('davidearn/burnout')", call. = FALSE)
}
if (!requireNamespace("gsl", quietly = TRUE)) {
  stop("Please install the gsl package: install.packages('gsl')", call. = FALSE)
}

get_script_dir <- function() {
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    if (!is.null(frames[[i]]$ofile)) {
      return(dirname(normalizePath(frames[[i]]$ofile, winslash = "/", mustWork = TRUE)))
    }
  }
  getwd()
}

script_dir <- get_script_dir()
main_file <- file.path(script_dir, "two_strain_burnout.R")
source(main_file)

output_dir <- file.path(script_dir, "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

## Parameters to scan.
R01_values <- seq(1.2, 3.0, by = 0.01)
R02_values <- seq(1.2, 3.0, by = 0.01)

## Parameters held fixed.
epsilon <- 0.01
N <- 1e6
k1 <- 1
k2 <- 1
dt <- 0.01
t_max <- 200
n_x <- 4000L
q_method <- "q_approx"

n_cores <- max(1L, parallel::detectCores() - 1L)

swap_entry_case <- function(x) {
  if (is.na(x)) return(NA_character_)
  x <- gsub("strain1", "STRAIN_ONE", x)
  x <- gsub("strain2", "strain1", x)
  gsub("STRAIN_ONE", "strain2", x)
}

make_row <- function(R01, R02, ans) {
  if (inherits(ans, "error")) {
    return(data.frame(
      R01 = R01, R02 = R02,
      boundary_event = "error",
      entry_case = NA_character_,
      Q1 = NA_real_, Q2 = NA_real_,
      surv1 = NA_real_, surv2 = NA_real_,
      both_burnout = NA_real_,
      strain1_only_burnout = NA_real_,
      strain2_only_burnout = NA_real_,
      neither_burnout = NA_real_,
      error_message = conditionMessage(ans)
    ))
  }

  if (!isTRUE(ans$boundary_event_found)) {
    return(data.frame(
      R01 = R01, R02 = R02,
      boundary_event = "no",
      entry_case = NA_character_,
      Q1 = NA_real_, Q2 = NA_real_,
      surv1 = NA_real_, surv2 = NA_real_,
      both_burnout = NA_real_,
      strain1_only_burnout = NA_real_,
      strain2_only_burnout = NA_real_,
      neither_burnout = NA_real_,
      error_message = NA_character_
    ))
  }

  Q1 <- unname(ans$strain_burnout["Q1"])
  Q2 <- unname(ans$strain_burnout["Q2"])
  out <- ans$outcomes_by_strain

  data.frame(
    R01 = R01, R02 = R02,
    boundary_event = ans$boundary_event,
    entry_case = ans$entry_case,
    Q1 = Q1,
    Q2 = Q2,
    surv1 = 1 - Q1,
    surv2 = 1 - Q2,
    both_burnout = unname(out["both_burnout"]),
    strain1_only_burnout = unname(out["strain1_only_burnout"]),
    strain2_only_burnout = unname(out["strain2_only_burnout"]),
    neither_burnout = unname(out["neither_burnout"]),
    error_message = NA_character_
  )
}

mirror_row <- function(row) {
  out <- row
  out$R01 <- row$R02
  out$R02 <- row$R01
  out$entry_case <- swap_entry_case(row$entry_case)
  out$Q1 <- row$Q2
  out$Q2 <- row$Q1
  out$surv1 <- row$surv2
  out$surv2 <- row$surv1
  out$strain1_only_burnout <- row$strain2_only_burnout
  out$strain2_only_burnout <- row$strain1_only_burnout
  out
}

run_one <- function(task) {
  i <- task[[1L]]
  j <- task[[2L]]
  R01 <- R01_values[i]
  R02 <- R02_values[j]

  ans <- tryCatch(
    two_strain_wave_burnout(
      R01 = R01,
      R02 = R02,
      epsilon = epsilon,
      N = N,
      k1 = k1,
      k2 = k2,
      dt = dt,
      t_max = t_max,
      n_x = n_x,
      q_method = q_method,
      store_trajectory = FALSE
    ),
    error = function(e) e
  )

  row <- make_row(R01, R02, ans)
  if (i == j) list(row) else list(row, mirror_row(row))
}

task_grid <- expand.grid(
  i = seq_along(R01_values),
  j = seq_along(R02_values),
  KEEP.OUT.ATTRS = FALSE
)
tasks <- split(task_grid[task_grid$i <= task_grid$j, ], seq_len(sum(task_grid$i <= task_grid$j)))

cat("Running ", length(tasks), " parameter points using ", n_cores, " core(s).\n", sep = "")

if (n_cores > 1L) {
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(
    cl,
    c(
      "main_file", "R01_values", "R02_values", "epsilon", "N", "k1", "k2",
      "dt", "t_max", "n_x", "q_method", "make_row", "mirror_row",
      "swap_entry_case", "run_one"
    ),
    envir = environment()
  )
  parallel::clusterEvalQ(cl, {
    if (!requireNamespace("burnout", quietly = TRUE)) stop("burnout package not found")
    if (!requireNamespace("gsl", quietly = TRUE)) stop("gsl package not found")
    source(main_file)
    NULL
  })
  row_chunks <- parallel::parLapply(cl, tasks, run_one)
} else {
  row_chunks <- lapply(tasks, run_one)
}

scan_results <- do.call(rbind, unlist(row_chunks, recursive = FALSE))
scan_results <- scan_results[order(scan_results$R01, scan_results$R02), ]
row.names(scan_results) <- NULL

csv_file <- file.path(output_dir, "two_strain_burnout_scan.csv")
write.csv(scan_results, csv_file, row.names = FALSE)

cat("\nSaved results to ", normalizePath(csv_file, winslash = "/", mustWork = FALSE), "\n", sep = "")