library(plagueMetapop)

## Endemic equilibrium for the deterministic ODE model.
##
## Force of infection: beta*I/K per susceptible per unit time.
## At endemic equilibrium (dI/dt = 0, I > 0): S* = gamma*K/beta = K/R0
## where R0 = beta/gamma.
##
## Jacobian analysis shows real eigenpart ~ -0.03/gen (logistic) and
## ~ -0.125/gen (linear), so by t = 500 gen the transient is negligible.
##
## Setup: one isolated patch (n_patch=1, alpha=0), strain 2 absent (I_init[2]=0),
## run to t = 500 disease generations.

beta  <- 2
gamma <- 1
K     <- 1e4
S_eq  <- gamma * K / beta   # = 5000; S* = K/R0 (R0 = beta/gamma = 2)

ode_base <- list(
  beta_vec  = c(beta, 0),
  K         = K,
  gamma     = c(gamma, gamma),
  r         = 0.125,
  n_patch   = 1,
  nt        = 500,
  alpha     = 0,
  I_init    = c(10, 0),
  def_file  = "ode_odin_def.R",
  stop_cond = NULL
)

## logistic vital dynamics (logistic_growth = 1, default)
run_logistic <- do.call(discrete_run, c(ode_base, list(logistic_growth = 1)))
S_final_logistic <- run_logistic |>
  dplyr::filter(step == max(step), state == "S") |>
  dplyr::pull(value)

expect_equal(S_final_logistic, S_eq, tolerance = 0.01,
             info = "ODE logistic demography: S converges to gamma*K/beta within 1%")

## linear restoring-force vital dynamics (logistic_growth = 0)
run_linear <- do.call(discrete_run,
                      modifyList(ode_base, list(logistic_growth = 0, r = 0.04)))
S_final_linear <- run_linear |>
  dplyr::filter(step == max(step), state == "S") |>
  dplyr::pull(value)

expect_equal(S_final_linear, S_eq, tolerance = 0.01,
             info = "ODE linear restoring-force demography: S converges to gamma*K/beta within 1%")
