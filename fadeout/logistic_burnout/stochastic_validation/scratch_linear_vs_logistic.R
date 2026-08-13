source(file.path("fadeout", "logistic_burnout", "logistic_burnout_functions.R"))

## --- Linear (Parsons et al.) model ---
## Classical no-demography SIR peak prevalence (eps -> 0 limit), Eq. 21:
##   ybar_0 = x_i - x_star*(1 + ln(x_i/x_star)), x_i = 1 (invasion into naive pop)
## Boundary layer: y_star = eps*(1 - 1/R0)
linear_ratio <- function(R0, eps) {
  x_star <- 1 / R0
  ybar0 <- 1 - x_star * (1 + log(1 / x_star))
  y_star <- eps * (1 - 1 / R0)
  list(ybar0 = ybar0, y_star = y_star, ratio = ybar0 / y_star)
}

cat("=== Linear (Parsons et al.) model: ybar_0 / y_star ===\n")
for (R0 in c(1.1, 1.2, 1.5)) {
  for (eps in c(0.1, 0.01, 0.001, 0.0007)) {
    res <- linear_ratio(R0, eps)
    cat(sprintf("R0=%.2f  eps=%-8s  ybar0=%.5f  y_star=%.6f  ratio=%8.2f\n",
                R0, eps, res$ybar0, res$y_star, res$ratio))
  }
  cat("\n")
}

## --- Our logistic model: same question, varying r instead of eps ---
cat("=== Logistic model: I_peak / I_BL, varying r ===\n")
for (R0 in c(1.1, 1.2, 1.5)) {
  for (r in c(0.1, 0.01, 0.001)) {
    K <- 100000  # large K so I_peak/I_BL is well resolved numerically; ratio is K-independent
    traj <- solve_logistic_SI(R0, r, K, I0 = 1, tmax = max(200, 20/r), dt = 0.01)
    entry <- find_logistic_boundary_entry(traj, R0, r, K)
    eq <- logistic_endemic_equilibrium(R0, r, K)
    if (isTRUE(entry$entry_found)) {
      ratio <- entry$I_peak / (K * eq$y_star)
      cat(sprintf("R0=%.2f  r=%-8s  y_peak=%.5f  y_star=%.6f  ratio=%8.2f\n",
                  R0, r, entry$I_peak / K, eq$y_star, ratio))
    } else {
      cat(sprintf("R0=%.2f  r=%-8s  status=%s\n", R0, r, entry$status))
    }
  }
  cat("\n")
}

## Classical (no-demography) SIR peak prevalence, for comparison with the r->0 limit:
cat("=== Classical SIR peak prevalence formula (r/eps -> 0 limit), for reference ===\n")
for (R0 in c(1.1, 1.2, 1.5)) {
  x_star <- 1 / R0
  ybar0 <- 1 - x_star * (1 + log(1 / x_star))
  cat(sprintf("R0=%.2f  ybar0(classical SIR) = %.5f\n", R0, ybar0))
}

## --- Does the ratio actually depend on K, fixing r and R0? ---
## y_star = r(R0-1)/R0^2 has no K in it at all. y_peak is a *density*
## computed from a trajectory that starts at I0=1 host, i.e. y(0)=1/K.
## So any K-dependence of y_peak/y_star must come in through this
## finite-population initial condition, not through the equilibrium formula.
cat("=== Fixed r, R0: does I_peak/I_BL (equivalently y_peak/y_star) depend on K? ===\n")
for (R0 in c(1.2, 1.5)) {
  for (r in c(0.05, 0.2)) {
    cat(sprintf("-- R0=%.2f  r=%.2f --\n", R0, r))
    for (K in c(1e3, 1e4, 1e5, 1e6, 1e7)) {
      traj <- solve_logistic_SI(R0, r, K, I0 = 1, tmax = max(200, 20 / r), dt = 0.01)
      entry <- find_logistic_boundary_entry(traj, R0, r, K)
      eq <- logistic_endemic_equilibrium(R0, r, K)
      if (isTRUE(entry$entry_found)) {
        ratio <- entry$I_peak / (K * eq$y_star)
        cat(sprintf("  K=%-9s y0=1/K=%-10.2e y_peak=%.6f  y_star=%.6f  ratio=%8.3f\n",
                    format(K, scientific = FALSE), 1 / K,
                    entry$I_peak / K, eq$y_star, ratio))
      } else {
        cat(sprintf("  K=%-9s status=%s\n", format(K, scientific = FALSE), entry$status))
      }
    }
  }
  cat("\n")
}

## --- Same finite-population check, but for the linear (Parsons et al.) model ---
## Normalized linear-vital-dynamics SI model:
##   dx/dt = eps*(1-x) - R0*x*y
##   dy/dt = (R0*x-1)*y
## Equilibrium: x_star=1/R0, y_star=eps*(1-1/R0) -- no N in the density formula,
## exactly parallel to the logistic y_star=r(R0-1)/R0^2 having no K in it.
## Any N-dependence of the peak density has to come from I0=1 => y(0)=1/N,
## just as I0=1 => y(0)=1/K did for the logistic model above.
solve_linear_SI <- function(R0, eps, N, I0 = 1, tmax = 50, dt = 0.02,
                             rtol = 1e-9, atol = 1e-11) {
  rhs <- function(time, state, parameters) {
    x <- state[["x"]]; y <- state[["y"]]
    list(c(
      x = eps * (1 - x) - R0 * x * y,
      y = (R0 * x - 1) * y
    ))
  }
  times <- sort(unique(c(seq(0, tmax, by = dt), tmax)))
  out <- as.data.frame(deSolve::ode(
    y = c(x = 1 - I0 / N, y = I0 / N), times = times,
    func = rhs, parms = NULL, method = "lsoda", rtol = rtol, atol = atol
  ))
  out$x <- pmax(out$x, 0); out$y <- pmax(out$y, 0)
  out
}

find_linear_peak <- function(trajectory, R0) {
  x_star <- 1 / R0
  idx <- which(head(trajectory$x, -1) > x_star & tail(trajectory$x, -1) <= x_star)
  if (!length(idx)) return(list(found = FALSE))
  i <- idx[1L]
  z0 <- trajectory$x[i]; z1 <- trajectory$x[i + 1L]
  frac <- (x_star - z0) / (z1 - z0)
  y_peak <- trajectory$y[i] + frac * (trajectory$y[i + 1L] - trajectory$y[i])
  list(found = TRUE, y_peak = y_peak)
}

cat("=== Linear model: fixed eps, R0, does y_peak/y_star depend on N? ===\n")
for (R0 in c(1.2, 1.5)) {
  for (eps in c(0.05, 0.2)) {
    y_star <- eps * (1 - 1 / R0)
    cat(sprintf("-- R0=%.2f  eps=%.2f  (y_star=%.6f) --\n", R0, eps, y_star))
    for (N in c(1e3, 1e4, 1e5, 1e6, 1e7)) {
      traj <- solve_linear_SI(R0, eps, N, I0 = 1, tmax = max(200, 20 / eps), dt = 0.01)
      pk <- find_linear_peak(traj, R0)
      if (isTRUE(pk$found)) {
        cat(sprintf("  N=%-9s y0=1/N=%-10.2e y_peak=%.6f  ratio=%8.3f\n",
                    format(N, scientific = FALSE), 1 / N, pk$y_peak, pk$y_peak / y_star))
      } else {
        cat(sprintf("  N=%-9s no peak found\n", format(N, scientific = FALSE)))
      }
    }
  }
  cat("\n")
}

## --- Does peak/boundary-layer ratio track where stochastic-validation error is large? ---
## Our validation grid fixed r=0.1 (see ../stochastic_validation/README.md). The
## documented large-error cells there are R0~1.1 (small K, multitrough truncation),
## R0=2.26426 (K=30000), and R0~3 (K=1000). Since the ratio is K-independent
## (established above), just scan R0 at r=0.1.
cat("=== Logistic model at our validation r=0.1: ratio across R0 ===\n")
r_fixed <- 0.1
K_probe <- 1e5
for (R0 in c(1.05, 1.1, 1.2, 1.5, 2, 2.26426, 2.5, 3, 4, 5)) {
  traj <- solve_logistic_SI(R0, r_fixed, K_probe, I0 = 1, tmax = max(200, 20 / r_fixed), dt = 0.01)
  entry <- find_logistic_boundary_entry(traj, R0, r_fixed, K_probe)
  eq <- logistic_endemic_equilibrium(R0, r_fixed, K_probe)
  if (isTRUE(entry$entry_found)) {
    ratio <- entry$I_peak / (K_probe * eq$y_star)
    cat(sprintf("R0=%-9s y_peak=%.5f  y_star=%.6f  ratio=%8.2f\n",
                R0, entry$I_peak / K_probe, eq$y_star, ratio))
  } else {
    cat(sprintf("R0=%-9s status=%s\n", R0, entry$status))
  }
}
cat("\n")

## --- Same scan for the linear model, at epsilon values from the actual
## Parsons et al. figure (talks/pix/parsons_burnout.png): x-axis "mean
## infectious period / mean lifetime" spans roughly 0 to 0.02 for n=1e6,
## with real diseases mostly clustered well below 0.001-0.01 and pneumonic
## plague itself essentially at eps -> 0. Our r=0.1 is far to the right of
## essentially their entire plotted range.
cat("=== Linear model at Parsons-et-al.-realistic epsilon: ratio across R0 ===\n")
for (eps in c(0.02, 0.005, 0.001, 0.0003, 0.0001)) {
  cat(sprintf("-- eps=%s --\n", eps))
  for (R0 in c(1.05, 1.1, 1.2, 1.5, 2, 2.26426, 2.5, 3, 4, 5)) {
    y_star <- eps * (1 - 1 / R0)
    ## The peak occurs on the epidemic (disease-generation) timescale, not the
    ## slow eps-driven recovery timescale, so cap tmax independent of eps.
    peak_tmax <- max(500, 60 / (R0 - 1))
    traj <- solve_linear_SI(R0, eps, N = 1e6, I0 = 1, tmax = peak_tmax, dt = 0.02)
    pk <- find_linear_peak(traj, R0)
    if (isTRUE(pk$found)) {
      cat(sprintf("  R0=%-9s y_peak=%.5f  y_star=%.6f  ratio=%10.2f\n",
                  R0, pk$y_peak, y_star, pk$y_peak / y_star))
    } else {
      cat(sprintf("  R0=%-9s no peak found\n", R0))
    }
  }
}
