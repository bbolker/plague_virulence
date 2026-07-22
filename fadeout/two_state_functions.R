## Helpers for the two-state occupancy closure with distinct infected loads in
## persistent and transient source patches.  These are deliberately separate
## from the original one-state approximation.

compute_transient_outbreak_summary <- function(R0, K, r, gamma = 1,
                                                I0 = 10, c = 0.5,
                                                integration_dt = 0.01) {
  values <- c(R0, K, r, gamma, I0, c, integration_dt)
  if (any(!is.finite(values)) || R0 <= 1 || K <= I0 || r <= 0 ||
      gamma <= 0 || I0 <= 0 || c <= 0 || integration_dt <= 0) {
    stop("Invalid deterministic transient-outbreak parameters")
  }

  ## This period closure is provisional.  The implemented population model has
  ## logistic susceptible recruitment rather than constant vital turnover, so
  ## mu=r and rho=gamma are an explicit approximation (as in the earlier
  ## two-state analysis), not an exact parameter identity.
  mu <- r
  rho <- gamma
  exact_argument <- mu * rho * (R0 - 1) - mu^2 * R0^2 / 4
  approximate_argument <- mu * rho * (R0 - 1)
  if (exact_argument <= 0) {
    return(list(
      oscillatory = FALSE, exact_argument = exact_argument,
      T_osc = NA_real_, T_osc_approx = 2 * pi / sqrt(approximate_argument),
      T_T = NA_real_, integral_I = NA_real_, Ibar_T = NA_real_,
      I_peak = NA_real_, time_of_I_peak = NA_real_, trajectory = NULL,
      mapping = "mu approximated by logistic growth parameter r; rho = gamma"
    ))
  }
  T_osc <- 2 * pi / sqrt(exact_argument)
  T_T <- c * T_osc
  times <- sort(unique(c(seq(0, T_T, by = integration_dt), T_T)))
  beta <- R0 * gamma

  ## These are the continuous-time equations corresponding to the repository's
  ## deterministic logistic single-patch model and plagueMetapop::ode_eq():
  ## dS/dt = r*S*(1-S/K) - beta*S*I/K; dI/dt = beta*S*I/K-gamma*I.
  rhs <- function(time, state, parameters) {
    S <- state[["S"]]
    I <- state[["I"]]
    incidence <- beta * S * I / K
    list(c(S = r * S * (1 - S / K) - incidence,
           I = incidence - gamma * I))
  }
  trajectory <- as.data.frame(deSolve::ode(
    y = c(S = K - I0, I = I0), times = times, func = rhs,
    parms = NULL, method = "lsoda", rtol = 1e-10, atol = 1e-12
  ))
  integral_I <- sum(diff(trajectory$time) *
    (head(trajectory$I, -1) + tail(trajectory$I, -1)) / 2)
  Ibar_T <- integral_I / T_T
  peak_index <- which.max(trajectory$I)
  I_peak <- trajectory$I[peak_index]
  time_of_I_peak <- trajectory$time[peak_index]
  if (!is.finite(Ibar_T) || Ibar_T <= 0 || Ibar_T > I_peak + 1e-7) {
    stop("Invalid deterministic transient infected-load summary")
  }
  list(
    oscillatory = TRUE, exact_argument = exact_argument,
    T_osc = T_osc, T_osc_approx = 2 * pi / sqrt(approximate_argument),
    T_T = T_T, integral_I = integral_I, Ibar_T = Ibar_T,
    I_peak = I_peak, time_of_I_peak = time_of_I_peak,
    trajectory = trajectory,
    mapping = "mu approximated by logistic growth parameter r; rho = gamma"
  )
}

solve_two_state_transient_Ibar <- function(time, p0, q0, P1, alpha,
                                           I_star, Ibar_T, T_T,
                                           numerical_tolerance = 1e-7) {
  inputs <- c(p0, q0, P1, alpha, I_star, Ibar_T, T_T)
  if (any(!is.finite(inputs)) || p0 < 0 || q0 < 0 || p0 + q0 > 1 ||
      P1 < 0 || P1 > 1 || alpha < 0 || I_star <= 0 || Ibar_T <= 0 ||
      T_T <= 0 || any(diff(time) < 0)) {
    stop("Invalid Ibar two-state parameters or initial conditions")
  }
  rhs <- function(t, state, parameters) {
    p <- state[["p"]]
    q <- state[["q"]]
    available <- max(0, 1 - p - q)
    B <- alpha * (I_star * p + Ibar_T * q) * available
    list(c(p = P1 * B, q = (1 - P1) * B - q / T_T))
  }
  ans <- as.data.frame(deSolve::ode(
    y = c(p = p0, q = q0), times = time, func = rhs, parms = NULL,
    method = "lsoda", rtol = 1e-9, atol = 1e-11
  ))
  ans$total <- ans$p + ans$q
  violation <- max(c(-ans$p, -ans$q, ans$total - 1), na.rm = TRUE)
  if (violation > numerical_tolerance) {
    stop("Substantive two-state state-space violation: ", violation)
  }
  ans$p[ans$p < 0] <- 0
  ans$q[ans$q < 0] <- 0
  ans$total <- ans$p + ans$q
  ans
}
