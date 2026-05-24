## modify basic vital-dynamics SIR equation:
##  * R0 (== beta) proportional to pop size
##  * set to discrete-time, single generation (all individuals recover after one time step)
##  * logistic growth (rounded) [no explicit death]
##  * for now, no recovered
##  * n_patch patches running in parallel

## based on https://mrc-ide.github.io/odin/articles/discrete.html

## Core equations for transitions between compartments:
update(S[]) <- S[i] - tot_incidence[i] + pop_change[i]
update(I[,1]) <- n_SI[i,1] + immig[i,1]
update(I[,2]) <- n_SI[i,2] + immig[i,2] + I2_seed[i]

## Individual probabilities of transition:
## mass action (scaled to carrying capacity)
hazard_SI[,] <- beta[j]*I[i,j]/K[i]
p_all[] <- -expm1(-sum(hazard_SI[i,]))
## probability that an infection is of type i
p_SI[,] <- if(sum(hazard_SI[i,])>0) hazard_SI[i,j]/sum(hazard_SI[i,]) else 1/n_strain

## Draws from binomial distributions for numbers changing between
## compartments:
## https://github.com/mrc-ide/sircovid/pull/178
## suggests first drawing a binomial for total number infected, then
##    multinomial from proportions
## odin rmultinom cannot accept array slices as arguments, so for
## n_strain=2 use sequential binomials (exact equivalent)
## FIXME: how much  more complicated to generalize this to n>2 strains?
##  (realistically, this is very low priority)
tot_incidence[] <- rbinom(S[i], p_all[i])
n_SI_strain1[]  <- rbinom(tot_incidence[i], p_SI[i,1])
n_SI[,1] <- n_SI_strain1[i]
n_SI[,2] <- tot_incidence[i] - n_SI_strain1[i]

## Vital dynamics (births/deaths): this is hokey
## (use constant births and logistic deaths instead?)
delta_log[] <- r[i]*S[i]*(1-S[i]/K[i])
pop_change[] <- if (delta_log[i] < 0) -rbinom(S[i], -delta_log[i]/S[i]) else rpois(delta_log[i])

## colonization
foi[] <- alpha*sum(I[,i])/n_patch
immig[,] <- rpois(foi[j])

## seed strain 2 at step == strain2_delay
I2_seed[] <- if (step == strain2_delay) I_ini[i,2] else 0
  
## Initial states:
initial(S[]) <- S_ini[i]
initial(I[,1]) <- I_ini[i,1]
initial(I[,2]) <- 0

## User defined parameters - default in parentheses:
n_patch <- user(100)
n_strain <- 2

## FIXME: add vector-valued defaults? (need to assign and then put into user(), functions not allowed)
beta[] <- user()
I_ini[,] <- user()
S_ini[] <- user()
alpha <- user(1e-5) ## between-patch transmission
strain2_delay <- user(0) ## steps before strain 2 is seeded
r[] <- user()   ## growth rate (per patch)
K[] <- user()   ## carrying capacity (per patch)

dim(S) <- n_patch
dim(I) <- c(n_patch, n_strain)
dim(hazard_SI) <- c(n_patch, n_strain)
dim(p_all) <- n_patch
dim(p_SI) <- c(n_patch, n_strain)
dim(tot_incidence) <- n_patch
dim(n_SI_strain1) <- n_patch
dim(n_SI) <- c(n_patch, n_strain)
dim(delta_log) <- n_patch
dim(pop_change) <- n_patch
dim(S_ini) <- n_patch
dim(I_ini) <- c(n_patch, n_strain)
dim(I2_seed) <- n_patch
dim(foi) <- n_strain
dim(immig) <- c(n_patch, n_strain)
dim(beta) <- n_strain
dim(r) <- n_patch
dim(K) <- n_patch
