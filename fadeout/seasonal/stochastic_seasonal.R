library(GillespieSSA2)
library(deSolve)

set.seed(12)

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
## Seasonal forcing
## ============================================================

## Seasonality is NOT sinusoidal.
##
## We use a piecewise-constant two-season model:
##
##   high-risk season: beta = beta_high
##   low-risk season:  beta = beta_low
##


season_period <- 365
high_duration <- 183

season_amp <- 0.2

beta_high <- beta * (1 + season_amp)
beta_low <- beta * (1 - season_amp)

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
cat("beta_high =", beta_high, "\n")
cat("beta_low =", beta_low, "\n")
cat("gamma =", gamma, "\n")
cat("mu =", mu, "\n")
cat("R0 =", R0, "\n")

## ============================================================
## Gillespie reactions
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

## Daily census output
dt <- 1

## ============================================================
## Seasonal boundaries
## ============================================================


boundaries <- sort(c(
  seq(
    from = high_duration,
    to = t_end,
    by = season_period
  ),
  seq(
    from = season_period,
    to = t_end,
    by = season_period
  )
))

segment_starts <- c(
  0,
  boundaries
)

segment_ends <- c(
  boundaries,
  t_end
)

keep <- segment_starts < segment_ends

segment_starts <- segment_starts[keep]
segment_ends <- segment_ends[keep]

## ============================================================
## Piecewise-constant seasonal Gillespie simulation
## ============================================================

current_state <- initial_state

stoch_time <- numeric(0)

stoch_state <- matrix(
  numeric(0),
  ncol = length(initial_state)
)

colnames(stoch_state) <- names(initial_state)

for (k in seq_along(segment_starts)) {
  
  segment_start <- segment_starts[k]
  segment_end <- segment_ends[k]
  
  ## Odd segments: high-risk season
  ## Even segments: low-risk season
  beta_current <- if (k %% 2 == 1) {
    beta_high
  } else {
    beta_low
  }
  
  segment_length <- segment_end - segment_start
  
  params_segment <- c(
    beta = beta_current,
    gamma = gamma,
    mu = mu,
    N = N
  )
  
  segment <- ssa(
    initial_state = current_state,
    reactions = reactions,
    params = params_segment,
    method = ssa_exact(),
    final_time = segment_length,
    census_interval = dt,
    verbose = FALSE
  )
  
  ## Convert local segment time to absolute time
  segment_time <- segment_start + segment$time
  
  ## Avoid duplicating seasonal boundary states
  if (length(stoch_time) > 0) {
    
    segment_time <- segment_time[-1]
    
    segment_state <- segment$state[
      -1,
      ,
      drop = FALSE
    ]
    
  } else {
    
    segment_state <- segment$state
  }
  
  ## Store results
  stoch_time <- c(
    stoch_time,
    segment_time
  )
  
  stoch_state <- rbind(
    stoch_state,
    segment_state
  )
  
  ## Final state becomes the initial state
  ## for the next seasonal segment
  current_state <- segment$state[
    nrow(segment$state),
    ,
    drop = TRUE
  ]
}

## ============================================================
## Deterministic seasonal demographic SIR
## ============================================================

sir_ode <- function(t, state, pars) {
  
  with(
    as.list(c(state, pars)),
    {
      
      ## Piecewise-constant seasonal transmission.
      ## This is NOT sinusoidal forcing.
      
      time_in_year <- t %% season_period
      
      beta_t <- if (time_in_year < high_duration) {
        beta_high
      } else {
        beta_low
      }
      
      dS <- mu * N -
        beta_t * S * I / N -
        mu * S
      
      dI <- beta_t * S * I / N -
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
  xlab = "Time (days)",
  ylab = "Infectives",
  main = "Deterministic seasonal SIR"
)

plot(
  stoch_time,
  stoch_state[, "I"],
  type = "l",
  lwd = 1,
  xlab = "Time (days)",
  ylab = "Infectives",
  main = "Stochastic seasonal SIR: Gillespie SSA"
)

par(mfrow = c(1, 1))

## ============================================================
## Post-transient comparison
## ============================================================

t_start <- 300

det_sub <- det[
  det$time >= t_start,
]

stoch_index <- stoch_time >= t_start

stoch_time_sub <- stoch_time[stoch_index]

stoch_I_sub <- stoch_state[
  stoch_index,
  "I"
]

ylim_post <- range(
  det_sub$I,
  stoch_I_sub
)

plot(
  det_sub$time,
  det_sub$I,
  type = "l",
  lwd = 1.5,
  ylim = ylim_post,
  xlab = "Time (days)",
  ylab = "Infectives",
  main = "Post-transient seasonal dynamics"
)

lines(
  stoch_time_sub,
  stoch_I_sub,
  lwd = 1
)

legend(
  "topright",
  legend = c(
    "Deterministic",
    "Stochastic"
  ),
  lty = 1,
  lwd = c(1.5, 1)
)