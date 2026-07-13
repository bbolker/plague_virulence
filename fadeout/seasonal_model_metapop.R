## Stochastic Euler metapopulation model
## with sinusoidal transmission seasonality
##
## Time unit: days
##
## beta[j] is the annual mean transmission rate.
##
## beta_eff[j](t) =
##   beta[j] *
##   (1 + seasonal_amp *
##        cos(2*pi*(t - peak_day)/season_period))
##
## Seasonal forcing is synchronous across patches.


## ------------------------------------------------------------
## State updates
## ------------------------------------------------------------

update(S[]) <- S[i] - tot_incidence[i] + pop_change[i]

update(I[,1]) <- I[i,1] + n_SI[i,1] - n_IR[i,1] + immig[i,1]

update(I[,2]) <- I[i,2] + n_SI[i,2] - n_IR[i,2] + immig[i,2] + I2_seed[i]


## ------------------------------------------------------------
## Seasonal transmission
## ------------------------------------------------------------

model_time <- step * dt

seasonal_multiplier <- 1 + seasonal_amp * cos(6.283185307179586 * (model_time - peak_day) / season_period)

beta_eff[] <- beta[i] * seasonal_multiplier


## ------------------------------------------------------------
## Infection hazards
## ------------------------------------------------------------

hazard_SI[,] <- beta_eff[j] * I[i,j] / K[i]

p_all[] <- -expm1(-sum(hazard_SI[i,]) * dt)

p_SI[,] <- if(sum(hazard_SI[i,]) > 0) hazard_SI[i,j] / sum(hazard_SI[i,]) else 1 / n_strain


## ------------------------------------------------------------
## Recovery / death
## ------------------------------------------------------------

p_recov[] <- -expm1(-gamma[i] * dt)


## ------------------------------------------------------------
## Infection draws
## ------------------------------------------------------------

tot_incidence[] <- rbinom(S[i], p_all[i])

n_SI_strain1[] <- rbinom(tot_incidence[i], p_SI[i,1])

n_SI[,1] <- n_SI_strain1[i]

n_SI[,2] <- tot_incidence[i] - n_SI_strain1[i]


## ------------------------------------------------------------
## Recovery / death draws
## ------------------------------------------------------------

n_IR[,] <- rbinom(I[i,j], p_recov[j])


## ------------------------------------------------------------
## Host demography
## ------------------------------------------------------------

delta_log[] <- r[i] * S[i] * (1 - S[i] / K[i]) * dt

pop_change[] <- if(delta_log[i] < 0) -rbinom(S[i], -delta_log[i] / S[i]) else rpois(delta_log[i])


## ------------------------------------------------------------
## Between-patch transmission
## ------------------------------------------------------------

foi[] <- alpha * sum(I[,i]) / n_patch

immig[,] <- rpois(foi[j] * dt)


## ------------------------------------------------------------
## Strain-2 delayed seed
## ------------------------------------------------------------

I2_active <- if(step == strain2_delay) 1.0 else 0.0

I2_seed[] <- I_ini[i,2] * I2_active


## ------------------------------------------------------------
## Initial states
## ------------------------------------------------------------

initial(S[]) <- S_ini[i]

initial(I[,1]) <- I_ini[i,1]

initial(I[,2]) <- I2_ini[i]


## ------------------------------------------------------------
## User parameters
## ------------------------------------------------------------

n_patch <- user(100)

n_strain <- 2

beta[] <- user()

gamma[] <- user()

dt <- user(1)

I_ini[,] <- user()

S_ini[] <- user()

I2_ini[] <- user()

alpha <- user(1e-4)

strain2_delay <- user(0L)

r[] <- user()

K[] <- user()

season_period <- user(365)

seasonal_amp <- user(0.3)

peak_day <- user(15)


## ------------------------------------------------------------
## Dimensions
## ------------------------------------------------------------

dim(S) <- n_patch

dim(I) <- c(n_patch, n_strain)

dim(beta) <- n_strain

dim(beta_eff) <- n_strain

dim(gamma) <- n_strain

dim(hazard_SI) <- c(n_patch, n_strain)

dim(p_all) <- n_patch

dim(p_SI) <- c(n_patch, n_strain)

dim(p_recov) <- n_strain

dim(tot_incidence) <- n_patch

dim(n_SI_strain1) <- n_patch

dim(n_SI) <- c(n_patch, n_strain)

dim(n_IR) <- c(n_patch, n_strain)

dim(delta_log) <- n_patch

dim(pop_change) <- n_patch

dim(S_ini) <- n_patch

dim(I_ini) <- c(n_patch, n_strain)

dim(I2_ini) <- n_patch

dim(I2_seed) <- n_patch

dim(foi) <- n_strain

dim(immig) <- c(n_patch, n_strain)

dim(r) <- n_patch

dim(K) <- n_patch
