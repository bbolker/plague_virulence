library(GillespieSSA2)
library(deSolve)

set.seed(123)

## ============================================================
## Parameters
## ============================================================

N <- 1e3

gamma <- 0.2
mu <- 0.02
R0 <- 2

beta <- R0 * (gamma + mu)

params <- c(
  beta = beta,
  gamma = gamma,
  mu = mu,
  N = N
)

## ============================================================
## Initial conditions
## ============================================================

I0 <- 10
R_init <- 0
S0 <- N - I0 - R_init

initial_state <- c(
  S = S0,
  I = I0,
  R = R_init
)

## ============================================================
## Print parameter information
## ============================================================

cat("N =", N, "\n")
cat("beta =", beta, "\n")
cat("gamma =", gamma, "\n")
cat("mu =", mu, "\n")
cat("R0 =", beta / (gamma + mu), "\n")

## ============================================================
## Gillespie stochastic demographic SIR
## ============================================================

reactions <- list(
  
  ## Birth into susceptible class
  reaction(
    propensity = ~ mu * N,
    effect = c(S = 1),
    name = "birth"
  ),
  
  ## Infection
  reaction(
    propensity = ~ beta * S * I / N,
    effect = c(S = -1, I = 1),
    name = "infection"
  ),
  
  ## Recovery
  reaction(
    propensity = ~ gamma * I,
    effect = c(I = -1, R = 1),
    name = "recovery"
  ),
  
  ## Susceptible death
  reaction(
    propensity = ~ mu * S,
    effect = c(S = -1),
    name = "S_death"
  ),
  
  ## Infective death
  reaction(
    propensity = ~ mu * I,
    effect = c(I = -1),
    name = "I_death"
  ),
  
  ## Recovered death
  reaction(
    propensity = ~ mu * R,
    effect = c(R = -1),
    name = "R_death"
  )
)

## ============================================================
## Simulation time
## ============================================================

t_end <- 1000
dt <- 0.5

## ============================================================
## Run Gillespie SSA
## ============================================================

stoch <- ssa(
  initial_state = initial_state,
  reactions = reactions,
  params = params,
  method = ssa_exact(),
  final_time = t_end,
  census_interval = dt,
  verbose = TRUE,
  sim_name = "Stochastic demographic SIR"
)

## ============================================================
## Deterministic demographic SIR
## ============================================================

sir_ode <- function(t, state, pars) {
  
  with(
    as.list(c(state, pars)),
    {
      
      dS <- mu * N -
        beta * S * I / N -
        mu * S
      
      dI <- beta * S * I / N -
        gamma * I -
        mu * I
      
      dR <- gamma * I -
        mu * R
      
      list(c(dS, dI, dR))
    }
  )
}

times <- seq(
  from = 0,
  to = t_end,
  by = dt
)

det <- ode(
  y = initial_state,
  times = times,
  func = sir_ode,
  parms = params
)

det <- as.data.frame(det)

## ============================================================
## Endemic equilibrium
## ============================================================

S_star <- N / R0

I_star <- mu * N * (R0 - 1) / beta

R_star <- N - S_star - I_star

cat("\nEndemic equilibrium:\n")
cat("S* =", S_star, "\n")
cat("I* =", I_star, "\n")
cat("R* =", R_star, "\n")

## ============================================================
## Plot infectives
## ============================================================

par(
  mfrow = c(2, 1),
  mar = c(4, 4, 3, 1)
)

plot(
  det$time,
  det$I,
  type = "l",
  lwd = 1.5,
  xlab = "Time",
  ylab = "Infectives",
  main = "Deterministic demographic SIR"
)

abline(
  h = I_star,
  lty = 2
)

plot(
  stoch$time,
  stoch$state[, "I"],
  type = "l",
  lwd = 1,
  xlab = "Time",
  ylab = "Infectives",
  main = "Stochastic demographic SIR: Gillespie SSA"
)

abline(
  h = I_star,
  lty = 2
)

par(mfrow = c(1, 1))

## ============================================================
## Plot post-transient infectives
## ============================================================

t_start <- 300

det_sub <- det[det$time >= t_start, ]

stoch_index <- stoch$time >= t_start
stoch_time <- stoch$time[stoch_index]
stoch_I <- stoch$state[stoch_index, "I"]

ylim_post <- range(
  det_sub$I,
  stoch_I
)

plot(
  det_sub$time,
  det_sub$I,
  type = "l",
  lwd = 1.5,
  ylim = ylim_post,
  xlab = "Time",
  ylab = "Infectives",
  main = "Post-transient dynamics"
)

lines(
  stoch_time,
  stoch_I,
  lwd = 1
)

abline(
  h = I_star,
  lty = 2
)

legend(
  "topright",
  legend = c("Deterministic", "Stochastic"),
  lty = 1,
  lwd = c(1.5, 1)
)