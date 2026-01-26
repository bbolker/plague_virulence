## based on https://mrc-ide.github.io/odin/articles/discrete.html

## Core equations for transitions between compartments:
## update(S) <- S - sum(n_SI[])
update(S) <- S - tot_incidence
update(I[]) <- I[i] + n_SI[i] - n_IR[i]
update(R) <- R + sum(n_IR)

## Individual probabilities of transition:
hazard_SI[] <- beta[i]*I[i]/N
p_all <- -expm1(-sum(hazard_SI))
p_SI[] <- if(sum(hazard_SI)>0) hazard_SI[i]/sum(hazard_SI) else 1/nstrains ##prevent NA when p_all=0
p_IR[] <- -expm1(-gamma[i]) ## 1 - exp(-gamma[i]) # I to R ## precompute?

## Draws from binomial distributions for numbers changing between
## compartments:
## https://github.com/mrc-ide/sircovid/pull/178
## suggests first drawing a binomial for total number infected, then
##    multinomial from proportions
## (seems easier than figuring out how to code a vector
##   [prob(no infection), prob(I1), prob(I2), ...])

tot_incidence <- rbinom(S, p_all)
n_SI[] <- rmultinom(tot_incidence, p_SI)
n_IR[] <- rbinom(I[i], p_IR[i])

## Initial states:
initial(S) <- N - sum(I_ini)
initial(I[]) <- I_ini[i]
initial(R) <- 0

## User defined parameters - default in parentheses:

nstrains <- 2
beta[] <- user()
gamma[] <- user()
I_ini[] <- user()
N <- user(10000)

dim(I) <- nstrains
dim(p_SI) <- nstrains
dim(hazard_SI) <- nstrains
dim(p_IR) <- nstrains
dim(n_IR) <- nstrains
dim(n_SI) <- nstrains
dim(I_ini) <- nstrains
dim(beta) <- nstrains
dim(gamma) <- nstrains
