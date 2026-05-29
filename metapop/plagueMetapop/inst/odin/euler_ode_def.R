## ODE (continuous-time) version of euler_odin_def.R using deriv():
##  * update(x) <- x + delta  replaced by  deriv(x) <- delta
##  * rpois(x) -> x;  rbinom(N, p) -> N*p  (deterministic means)
##  * p_all and p_recov are now instantaneous rates per unit time,
##    not per-step probabilities (expm1 removed; no explicit dt)
##  * reedfrost and strain2_delay are not applicable to ODE models and are dropped;
##    strain 2 starts at t=0 from I_ini[,2] (set to 0 for "strain 2 absent initially")

## Core derivatives:
deriv(S[])   <- -tot_incidence[i] + pop_change[i]
deriv(I[,1]) <- n_SI[i,1] - n_IR[i,1] + immig[i,1]
deriv(I[,2]) <- n_SI[i,2] - n_IR[i,2] + immig[i,2]

## Force of infection (instantaneous rate per unit time):
hazard_SI[,] <- beta[j]*I[i,j]/K[i]
p_all[]      <- sum(hazard_SI[i,])
## probability that a new infection is of strain j (ratio of rates; unchanged)
p_SI[,] <- if(sum(hazard_SI[i,])>0) hazard_SI[i,j]/sum(hazard_SI[i,]) else 1/n_strain

## Recovery rate per unit time:
p_recov[] <- gamma[i]

## New infections (deterministic flows):
tot_incidence[]  <- p_all[i] * S[i]
n_SI_strain1[]   <- tot_incidence[i] * p_SI[i,1]
n_SI[,1] <- n_SI_strain1[i]
n_SI[,2] <- tot_incidence[i] - n_SI_strain1[i]

## Recovery (deterministic flow):
n_IR[,] <- p_recov[j] * I[i,j]

## Vital dynamics (rate per unit time; no dt):
delta_log[] <- if (logistic_growth == 1) r[i]*S[i]*(1-S[i]/K[i]) else r[i]*(K[i]-S[i])
pop_change[] <- delta_log[i]

## Immigration (rate per unit time; no dt):
foi[]    <- alpha * sum(I[,i]) / n_patch
immig[,] <- foi[j]

## Initial states:
initial(S[])   <- S_ini[i]
initial(I[,1]) <- I_ini[i,1]
initial(I[,2]) <- I_ini[i,2]

## User-defined parameters:
n_patch  <- user(100)
n_strain <- 2

beta[]   <- user()
gamma[]  <- user()
I_ini[,] <- user()
S_ini[]  <- user()
alpha           <- user(1e-5)
r[]             <- user()
K[]             <- user()
logistic_growth <- user(1)  ## 1 = standard logistic; 0 = linear restoring force

dim(S)             <- n_patch
dim(I)             <- c(n_patch, n_strain)
dim(hazard_SI)     <- c(n_patch, n_strain)
dim(p_all)         <- n_patch
dim(p_SI)          <- c(n_patch, n_strain)
dim(p_recov)       <- n_strain
dim(tot_incidence) <- n_patch
dim(n_SI_strain1)  <- n_patch
dim(n_SI)          <- c(n_patch, n_strain)
dim(n_IR)          <- c(n_patch, n_strain)
dim(delta_log)     <- n_patch
dim(pop_change)    <- n_patch
dim(S_ini)         <- n_patch
dim(I_ini)         <- c(n_patch, n_strain)
dim(foi)           <- n_strain
dim(immig)         <- c(n_patch, n_strain)
dim(beta)          <- n_strain
dim(gamma)         <- n_strain
dim(r)             <- n_patch
dim(K)             <- n_patch
