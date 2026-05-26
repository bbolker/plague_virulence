## Deterministic version of euler_odin_def.R:
##  * rpois(lambda) replaced by lambda
##  * rbinom(N, p)  replaced by N*p
## State variables become real-valued (not integer).

## Core equations for transitions between compartments:
update(S[]) <- S[i] - tot_incidence[i] ## + pop_change[i]
update(I[,1]) <- I[i,1] + n_SI[i,1] - n_IR[i,1] + immig[i,1]
update(I[,2]) <- I[i,2] + n_SI[i,2] - n_IR[i,2] + immig[i,2] + I2_seed[i]

## Individual probabilities of transition:
## mass action (scaled to carrying capacity), multiplied by dt
hazard_SI[,] <- beta[j]*I[i,j]/K[i]
p_all[] <- -expm1(-sum(hazard_SI[i,])*dt)
## probability that a new infection is of strain j
p_SI[,] <- if(sum(hazard_SI[i,])>0) hazard_SI[i,j]/sum(hazard_SI[i,]) else 1/n_strain

## Recovery/death probability per strain (patch-independent)
p_recov[] <- -expm1(-gamma[i]*dt)

## Deterministic flow for new infections:
tot_incidence[]  <- S[i] * p_all[i]
n_SI_strain1[]   <- tot_incidence[i] * p_SI[i,1]
n_SI[,1] <- n_SI_strain1[i]
n_SI[,2] <- tot_incidence[i] - n_SI_strain1[i]

## Deterministic recovery/death flow
n_IR[,] <- I[i,j] * p_recov[j]

## Vital dynamics (births/deaths), scaled by dt
delta_log[] <- r[i]*S[i]*(1-S[i]/K[i])*dt
pop_change[] <- delta_log[i]

## colonization (rate scaled by dt)
foi[] <- alpha*sum(I[,i])/n_patch
immig[,] <- foi[j]*dt

## seed strain 2 at step == strain2_delay
## separate scalar indicator from array equation to avoid odin if/step/array interaction
I2_active <- if (step == strain2_delay) 1.0 else 0.0
I2_seed[] <- I_ini[i,2] * I2_active

## Initial states:
initial(S[]) <- S_ini[i]
initial(I[,1]) <- I_ini[i,1]
initial(I[,2]) <- 0

## User defined parameters - default in parentheses:
n_patch <- user(100)
n_strain <- 2

beta[] <- user()
gamma[] <- user()
dt <- user(1)
I_ini[,] <- user()
S_ini[] <- user()
alpha <- user(1e-5) ## between-patch transmission
strain2_delay <- user(0L) ## steps before strain 2 is seeded (integer: compared to step)
r[] <- user()   ## growth rate (per patch)
K[] <- user()   ## carrying capacity (per patch)

dim(S) <- n_patch
dim(I) <- c(n_patch, n_strain)
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
dim(I2_seed) <- n_patch
dim(foi) <- n_strain
dim(immig) <- c(n_patch, n_strain)
dim(beta) <- n_strain
dim(gamma) <- n_strain
dim(r) <- n_patch
dim(K) <- n_patch
