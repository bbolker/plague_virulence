## Scratch: exact K -> Inf limit of y_peak via unstable-manifold seeding.
##
## The (x,y) ODE dx/dt=rx(1-x)-R0xy, dy/dt=(R0x-1)y contains no K at all;
## K only enters through the initial condition x(0)=1-I0/K, y(0)=I0/K. As
## K -> Inf this initial point -> (1,0), which is itself a fixed point (the
## disease-free state) -- a saddle, with Jacobian
##   J(1,0) = [[-r, -R0], [0, R0-1]]
## eigenvalues -r (stable) and R0-1 (unstable). Any I0=1-with-finite-K
## trajectory has a nonzero component along the unstable eigenvector, so for
## large K it tracks the unique unstable-manifold trajectory emanating from
## (1,0) (the exponentially-growing component dominates the exponentially-
## decaying one). So instead of brute-force large K, seed the ODE a tiny
## distance from (1,0) along the exact unstable eigenvector direction --
## this converges to the true K=Inf limit much faster and more precisely.
source(file.path("fadeout", "logistic_burnout", "logistic_burnout_functions.R"))

unstable_eigvec <- function(R0, r) {
  ## J(1,0) = [[-r,-R0],[0,R0-1]]; eigenvector for eigenvalue (R0-1):
  ## (-r-(R0-1)) vx - R0 vy = 0  =>  vy = ((1-r-R0)/R0) vx
  vx <- -1
  vy <- ((1 - r - R0) / R0) * vx
  v <- c(vx, vy) / sqrt(vx^2 + vy^2)
  v
}

Kinf_peak <- function(R0, r, delta = 1e-9, tmax = 400, dt = 0.01) {
  v <- unstable_eigvec(R0, r)
  x0 <- 1 + delta * v[1]
  y0 <- delta * v[2]
  rhs <- function(time, state, parameters) {
    x <- state[["x"]]; y <- state[["y"]]
    list(c(x = r * x * (1 - x) - R0 * x * y, y = (R0 * x - 1) * y))
  }
  times <- sort(unique(c(seq(0, tmax, by = dt), tmax)))
  traj <- as.data.frame(deSolve::ode(
    y = c(x = x0, y = y0), times = times, func = rhs, parms = NULL,
    method = "lsoda", rtol = 1e-11, atol = 1e-13
  ))
  xstar <- 1 / R0
  idx <- which(head(traj$x, -1) > xstar & tail(traj$x, -1) <= xstar)
  if (!length(idx)) return(list(found = FALSE))
  i <- idx[1L]
  frac <- (xstar - traj$x[i]) / (traj$x[i + 1L] - traj$x[i])
  y_peak <- traj$y[i] + frac * (traj$y[i + 1L] - traj$y[i])
  list(found = TRUE, y_peak = y_peak)
}

cat("=== K=Inf limit (unstable-manifold seeding) vs finite-K numerics ===\n")
for (R0 in c(1.1, 1.2, 1.5, 2, 3, 5)) {
  r <- 0.1
  lim <- Kinf_peak(R0, r)
  ## compare against explicit finite-K solves
  finite <- sapply(c(1e3, 1e4, 1e5, 1e6, 1e7), function(K) {
    traj <- solve_logistic_SI(R0, r, K, I0 = 1, tmax = 200, dt = 0.01)
    entry <- find_logistic_boundary_entry(traj, R0, r, K)
    if (isTRUE(entry$entry_found)) entry$I_peak / K else NA_real_
  })
  cat(sprintf("R0=%.1f: K=Inf limit y_peak=%.6f | K=1e3:%.6f 1e4:%.6f 1e5:%.6f 1e6:%.6f 1e7:%.6f\n",
              R0, lim$y_peak, finite[1], finite[2], finite[3], finite[4], finite[5]))
}

## Also check sensitivity of the K=Inf estimate itself to the seed distance delta
cat("\n=== K=Inf estimate sensitivity to seed distance delta (R0=1.2, r=0.1) ===\n")
for (delta in c(1e-3, 1e-5, 1e-7, 1e-9, 1e-11)) {
  lim <- Kinf_peak(1.2, 0.1, delta = delta)
  cat(sprintf("delta=%s -> y_peak=%.8f\n", delta, lim$y_peak))
}
