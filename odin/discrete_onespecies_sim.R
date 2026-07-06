library(plagueMetapop)
library(future)
library(ggplot2); theme_set(theme_bw())

nsim <- 20
ncores <- 10

dd <- expand.grid(
    R0 = seq(1.1, 3.0, by = 0.1),
    alpha = 10^seq(-5.5, -3.5, by = 0.5),  
    K = 10^seq(4:7)
)
nrow(dd) ## 400

plan(multicore(workers = 10))
## plan(sequential)
set.seed(101)

ctr <- 0
nsim <- 5
FUN <- function(x) {
    out <- with(x,
         discrete_run(alpha = alpha,
                      K = K,
                      beta_vec = c(R0, 1),
                      I_init = c(10, 0),
                      nt = 500))
    sumfun_discrete(out)
}
FUN(dd[1,])

    
apply(dd, 1,
