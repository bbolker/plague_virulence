## Scratch: peak/boundary-layer ratio vs R0, logistic (r=0.1, our validation
## grid) overlaid with linear/Parsons-et-al. model at several epsilon values
## (0.1 matched to r=0.1 for direct comparison, plus realistic small-eps
## values read off talks/pix/parsons_burnout.png). Companion to
## scratch_linear_vs_logistic.R and scratch_peak_boundary_ratio_vs_R0.R.
## Not part of the committed pipeline.
source(file.path("fadeout", "logistic_burnout", "logistic_burnout_functions.R"))

## --- linear (Parsons et al.) deterministic model, matching
## scratch_linear_vs_logistic.R's definitions ---
solve_linear_SI <- function(R0, eps, N, I0 = 1, tmax = 50, dt = 0.02,
                             rtol = 1e-9, atol = 1e-11) {
  rhs <- function(time, state, parameters) {
    x <- state[["x"]]; y <- state[["y"]]
    list(c(x = eps * (1 - x) - R0 * x * y, y = (R0 * x - 1) * y))
  }
  times <- sort(unique(c(seq(0, tmax, by = dt), tmax)))
  out <- as.data.frame(deSolve::ode(
    y = c(x = 1 - I0 / N, y = I0 / N), times = times,
    func = rhs, parms = NULL, method = "lsoda", rtol = rtol, atol = atol
  ))
  out$x <- pmax(out$x, 0); out$y <- pmax(out$y, 0)
  out
}

find_linear_ratio <- function(R0, eps, N = 1e6, initial_tmax = 200,
                               maximum_tmax = 12800, dt = 0.02) {
  xstar <- 1 / R0
  y_star <- eps * (1 - 1 / R0)
  horizon <- initial_tmax
  repeat {
    traj <- solve_linear_SI(R0, eps, N, I0 = 1, tmax = horizon, dt = dt)
    idx <- which(head(traj$x, -1) > xstar & tail(traj$x, -1) <= xstar)
    if (length(idx) || horizon >= maximum_tmax) break
    horizon <- min(maximum_tmax, 2 * horizon)
  }
  if (!length(idx)) {
    return(data.frame(R0 = R0, eps = eps, y_peak = NA_real_, y_BL = y_star,
                       ratio = NA_real_, status = "no_epidemic_peak", horizon_used = horizon))
  }
  i <- idx[1L]
  frac <- (xstar - traj$x[i]) / (traj$x[i + 1L] - traj$x[i])
  y_peak <- traj$y[i] + frac * (traj$y[i + 1L] - traj$y[i])
  data.frame(R0 = R0, eps = eps, y_peak = y_peak, y_BL = y_star,
             ratio = y_peak / y_star, status = "success", horizon_used = horizon)
}

find_logistic_ratio <- function(R0, r, K = 1e5, initial_tmax = 100,
                                 maximum_tmax = 6400, dt = 0.02) {
  eq <- logistic_endemic_equilibrium(R0, r, K)
  horizon <- initial_tmax
  repeat {
    traj <- solve_logistic_SI(R0, r, K, I0 = 1, tmax = horizon, dt = dt)
    entry <- find_logistic_boundary_entry(traj, R0, r, K)
    if (entry$entry_found || horizon >= maximum_tmax) break
    horizon <- min(maximum_tmax, 2 * horizon)
  }
  if (!isTRUE(entry$entry_found)) {
    return(data.frame(R0 = R0, r = r, y_peak = NA_real_, y_BL = eq$y_star,
                       ratio = NA_real_, status = entry$status, horizon_used = horizon))
  }
  data.frame(R0 = R0, r = r, y_peak = entry$I_peak / K, y_BL = eq$y_star,
             ratio = (entry$I_peak / K) / eq$y_star, status = "success", horizon_used = horizon)
}

## Common R0 grid: log-spaced in (R0-1) from near 1 to 4, i.e. R0 in (1, 5].
R0_grid <- 1 + exp(seq(log(0.003), log(4), length.out = 26))

out_dir <- file.path("fadeout", "logistic_burnout", "stochastic_validation", "outputs")

logistic_res <- do.call(rbind, lapply(R0_grid, find_logistic_ratio, r = 0.1))
logistic_res$model <- "logistic"
logistic_res$param <- "r=0.1"
write.csv(logistic_res, file.path(out_dir, "scratch_logistic_ratio_full.csv"), row.names = FALSE)

eps_values <- c(0.1, 0.02, 0.001, 0.0001)
linear_list <- list()
for (eps in eps_values) {
  cat(sprintf("linear model, eps=%s\n", eps))
  res <- do.call(rbind, lapply(R0_grid, find_linear_ratio, eps = eps))
  res$model <- "linear"
  res$param <- sprintf("eps=%s", eps)
  linear_list[[length(linear_list) + 1L]] <- res
}
linear_res <- do.call(rbind, linear_list)
write.csv(linear_res, file.path(out_dir, "scratch_linear_ratio_full.csv"), row.names = FALSE)

cat("\nDone. Rows with status=success per curve:\n")
cat("logistic r=0.1: ", sum(logistic_res$status == "success"), "/", nrow(logistic_res), "\n")
for (eps in eps_values) {
  sub <- linear_res[linear_res$eps == eps, ]
  cat(sprintf("linear eps=%s: %d/%d\n", eps, sum(sub$status == "success"), nrow(sub)))
}
