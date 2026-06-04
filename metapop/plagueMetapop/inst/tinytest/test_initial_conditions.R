library(plagueMetapop)

## Verify that the first row of run_simulator_odin output (step/t == 0)
## matches the initial conditions supplied to make_simulator_odin.
## n_patch = 1 gives deterministic initial conditions (no Poisson draws).

K      <- 1e4
I_init <- c(20, 15)

## S_ini = K - round(I_init[1]) - round(I_init[2]) for all models
## (make_simulator_odin computes S_ini = K - rowSums(I_ini_mat))
I1_exp <- round(I_init[1])
S_exp  <- K - round(I_init[1]) - round(I_init[2])

## --- discrete model ---
## I[,2] initial value comes from I2_ini = 0, NOT from I_init[2].
## Strain 2 is seeded later via strain2_delay / I2_seed[].
mod_disc <- make_simulator_odin(
  beta_vec = c(1.5, 2.5),
  K        = K,
  n_patch  = 1L,
  nt       = 10L,
  I_init   = I_init,
  alpha    = 0
)
raw_disc <- run_simulator_odin(mod_disc, stop_cond = NULL)
row0     <- raw_disc[raw_disc[, "step"] == 0, , drop = FALSE]

expect_equal(as.numeric(row0[, "S[1]"]),   S_exp,  info = "discrete: S[1] at step 0 matches S_ini")
expect_equal(as.numeric(row0[, "I[1,1]"]), I1_exp, info = "discrete: I[1,1] at step 0 matches I_ini[1,1]")
expect_equal(as.numeric(row0[, "I[1,2]"]), 0,      info = "discrete: I[1,2] at step 0 is 0 (I2_ini; seeded later)")

## --- ODE model ---
## ODE uses initial(I[,2]) <- I_ini[i,2] directly, so strain 2 starts at
## round(I_init[2]) rather than 0.
I2_exp_ode <- round(I_init[2])

mod_ode <- make_simulator_odin(
  beta_vec = c(1.5, 2.5),
  K        = K,
  gamma    = c(1, 1),
  r        = 0.125,
  n_patch  = 1L,
  nt       = 10L,
  I_init   = I_init,
  alpha    = 0,
  def_file = "ode_odin_def.R"
)
raw_ode <- run_simulator_odin(mod_ode, stop_cond = NULL)
row_t0  <- raw_ode[1, , drop = FALSE]   ## first row is always t=0

expect_equal(as.numeric(row_t0[, "S[1]"]),   S_exp,      info = "ODE: S[1] at t=0 matches S_ini")
expect_equal(as.numeric(row_t0[, "I[1,1]"]), I1_exp,     info = "ODE: I[1,1] at t=0 matches I_ini[1,1]")
expect_equal(as.numeric(row_t0[, "I[1,2]"]), I2_exp_ode, info = "ODE: I[1,2] at t=0 matches I_ini[1,2] (not 0)")
