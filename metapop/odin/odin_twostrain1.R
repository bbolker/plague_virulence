## modify basic vital-dynamics SIR equation:
##  * R0 (== beta) proportional to pop size
##  * set to discrete-time, single generation (all individuals recover after one time step)
##  * logistic growth (rounded) [no explicit death]
##  * for now, no recovered

## based on https://mrc-ide.github.io/odin/articles/discrete.html

## Core equations for transitions between compartments:
## update(S) <- S - sum(n_SI[])

update(S) <- S - tot_incidence + n_births - n_deaths_S
update(I[]) <- n_SI[i]

## Individual probabilities of transition:
hazard_SI[] <- beta[i]*I[i]/N
p_all <- -expm1(-sum(hazard_SI))
p_SI[] <- if(sum(hazard_SI)>0) hazard_SI[i]/sum(hazard_SI) else 1/nstrains ##prevent NA when p_all=0

## Draws from binomial distributions for numbers changing between
## compartments:
## https://github.com/mrc-ide/sircovid/pull/178
## suggests first drawing a binomial for total number infected, then
##    multinomial from proportions
## (seems easier than figuring out how to code a vector
##   [prob(no infection), prob(I1), prob(I2), ...])

tot_incidence <- rbinom(S, p_all)
n_SI[] <- rmultinom(tot_incidence, p_SI)

## Vital dynamics (births/deaths): this is hokey
## (use constant births and logistic deaths instead?)
delta_log <- r*N*(1-N/K)
if (delta_log < 0) {
  pop_change <- -rbinom(S, -delta_log/S)
} else {
  pop_change <- rpois(delta_log)
}

## Initial states:
initial(S) <- S_ini
initial(I[]) <- I_ini[i]
initial(R) <- N-S_ini-sum(I_ini)

## User defined parameters - default in parentheses:

nstrains <- 2
beta[] <- user()
gamma[] <- user()
I_ini[] <- user()
S_ini <- user()
N <- user(10000)
mu <- user(0.0)  ## birth/death rate 

dim(I) <- nstrains
dim(p_SI) <- nstrains
dim(hazard_SI) <- nstrains
dim(p_IR) <- nstrains
dim(n_IR) <- nstrains
dim(n_SI) <- nstrains
dim(I_ini) <- nstrains
dim(beta) <- nstrains
dim(gamma) <- nstrains
dim(n_deaths_I) <- nstrains
