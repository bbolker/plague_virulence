## Scratch: peak/boundary-layer ratio y_peak/y_BL as a function of R0 at the
## validation grid's fixed r=0.1, extended down toward the overdamped
## threshold (below which no deterministic epidemic peak exists at all).
## Companion to scratch_linear_vs_logistic.R. Not part of the committed
## logistic_burnout pipeline.
source(file.path("fadeout", "logistic_burnout", "logistic_burnout_functions.R"))

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")

r_fixed <- 0.1

## Overdamped/critical-damping threshold: x(t) fails to cross x*=1/R0 (and so
## no peak exists) when r(R0-1) < r^2 R0^2/4, i.e. R0 < R0_crit(r), where
## R0_crit solves r*R0^2 - 4*R0 + 4 = 0:
R0_crit <- function(r) 2 * (1 - sqrt(1 - r)) / r
cat(sprintf("R0_crit(r=%.3f) = %.6f  (documented in ../stochastic_validation/README.md as the overdamped boundary)\n",
            r_fixed, R0_crit(r_fixed)))

find_ratio <- function(R0, r, K = 1e5, initial_tmax = 100, maximum_tmax = 6400, dt = 0.02) {
  ## NOTE: unlike logistic_burnout_probability()'s retry loop (which only
  ## re-extends the horizon when a peak WAS found but no post-peak boundary
  ## crossing was -- a documented limitation in ../stochastic_validation/README.md),
  ## this one also re-extends on "no_epidemic_peak", since near R0_crit the
  ## approach to x* can be slow enough that the peak simply hasn't happened
  ## yet within a short horizon (as opposed to genuinely never happening,
  ## which only occurs for R0 < R0_crit(r)).
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

## Dense grid pushing down toward R0_crit, plus the existing validation-grid
## reference points (including the documented 2.26426 log-grid point).
crit <- R0_crit(r_fixed)
R0_near_crit <- crit + c(0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.03, 0.05)
R0_reference <- c(1.05, 1.1, 1.2, 1.5, 2, 2.26426, 2.5, 3, 4, 5)
R0_grid <- sort(unique(c(R0_near_crit, R0_reference)))

results <- do.call(rbind, lapply(R0_grid, find_ratio, r = r_fixed))
print(results, row.names = FALSE)

out_dir <- file.path("fadeout", "logistic_burnout", "stochastic_validation", "outputs")
write.csv(results, file.path(out_dir, "scratch_peak_boundary_ratio_vs_R0.csv"), row.names = FALSE)

## Also show r-sensitivity at a couple of R0 near the low end, to see how
## R0_crit itself moves.
cat("\nR0_crit at other r values used in the stochastic_validation r_sensitivity grid:\n")
for (r in c(0.05, 0.1, 0.125, 0.2)) {
  cat(sprintf("  r=%.3f -> R0_crit=%.5f\n", r, R0_crit(r)))
}
