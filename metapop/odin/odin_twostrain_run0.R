library(odin)
library(dde) ## odin insists on this

nt <- 1000
beta_vec <- c(1.5, 2.5)
r <- 0.125
K <- 1e6
I_init <- c(10, 10)

odin_fn <- here::here("metapop/odin", "odin_twostrain1.R")
twostrain_generator <- odin::odin(odin_fn)

odin_fn <- twostrain_generator$new(beta = c(beta_vec[1], beta_vec[2]),
                                   r = r,
                                   K = K,
                                   I_ini = I_init,
                                   S_ini = K - sum(I_init))
res <- odin_fn$run(0:nt)


