## Single-strain Euler SIR model (n_patch patches).
## Simplified from euler_odin_def.R: all strain-2 machinery removed,
## beta and gamma are scalars, I_ini is a 1-D array.

## State updates
update(S[]) <- S[i] - n_SI[i] + pop_change[i]
update(I[]) <- I[i] + n_SI[i] - n_IR[i] + immig[i]

## Infection hazard and probability (hazard scaled by dt)
hazard_SI[] <- beta * I[i] / K[i]
p_SI[]      <- -expm1(-hazard_SI[i] * dt)

## Recovery probability (patch-independent)
p_recov <- -expm1(-gamma * dt)

## Stochastic draws
n_SI[] <- rbinom(S[i], p_SI[i])
n_IR[] <- rbinom(I[i], p_recov)

## Vital dynamics (logistic, scaled by dt)
delta_log[]  <- r[i] * S[i] * (1 - S[i] / K[i]) * dt
pop_change[] <- if (delta_log[i] < 0) -rbinom(S[i], -delta_log[i] / S[i]) else rpois(delta_log[i])

## Between-patch colonization (rate scaled by dt)
foi     <- alpha * sum(I[]) / n_patch
immig[] <- rpois(foi * dt)

## Initial conditions
initial(S[]) <- S_ini[i]
initial(I[]) <- I_ini[i]

## Parameters
n_patch <- user(100)

beta    <- user()
gamma   <- user(1)
dt      <- user(1)
I_ini[] <- user()
S_ini[] <- user()
alpha   <- user(1e-5)
r[]     <- user()
K[]     <- user()

dim(S)         <- n_patch
dim(I)         <- n_patch
dim(hazard_SI) <- n_patch
dim(p_SI)      <- n_patch
dim(n_SI)      <- n_patch
dim(n_IR)      <- n_patch
dim(delta_log) <- n_patch
dim(pop_change) <- n_patch
dim(S_ini)     <- n_patch
dim(I_ini)     <- n_patch
dim(immig)     <- n_patch
dim(r)         <- n_patch
dim(K)         <- n_patch
