## avoid T, it's a built-in synonym for TRUE
#' @param tt number of time steps
#' @param B0 infection rate between patches
#' @param r  growth rate of rat population
#' @param D  death rate due to pathogen
#' @param epsilon ratio of pathogen to host generation time
#' @param R0 intrinsic reproductive number of pathogen
#' @param ki ?? for infected patches
#' @param ks ?? for susc patches
#' @param c0
#' @param K  rat carrying capacity per patch
#' @param N number of patches (only affects discrete-extinction logic)
#' @param starting conditions
simfun <- function(tt  = 100,
                   B0 = 0.5,
                   r = 0.5,
                   D = 0.5,
                   epsilon=0.05,
                   R0=10,
                   c0=0.2,
                   K=300,
                   S = K,       ##fixed susceptible patch population
                   
                   start = c(p=0.01, N=K)) {
  
  z=final_size(R0)
  N_DFE=c0*(1+r)*(1-z*D)*S/(1-(1+r)*(1-z*D)*(1-c0)) 
  cat("NI at DFE is",N_DFE,"\n")
  
  if (!require("burnout")) stop("please install the 'burnout' package: ",
                                "`remotes::install_github('davidearn/burnout')`")
  
  ## Initialize parameters
  P <- N <- rep(NA, tt)  ## P for fraction of infected patches, N for population density of infected patches
  
  P[1] <- start[["p"]]
  N[1] <- start[["N"]]
  
  ## Simulation loop
  for (t in 1:(tt-1)) {
    ## Current state
    p = P[t]
    n = N[t]
    
    ##ignore stocastic extinction
    
    ## host movement and infection between patches(caused by host movement)
    p = p + B0*p*(1-p)                 ##infection of S patches (constant B for simplicity)
    n <- n + c0*(S-n)*(1-p)            ## host (rat) movement between patches
    ## ignore ni change of n because S patches convert to I (combined with movement event for simplicity)
      
    ##epidemic:death of hosts and burnout
    ##death of hosts  
    n <- n * (1 - z*D)  
    ##burnout
    P1=P1_prob(R0,epsilon,k=1,N=n)
    p <-P1*p 
    
    ## logistic growth
    #n=n+r*n*(1-n/K)
    
    ## it doesn't have to be logistic if death proportion is high(?)
    n=n*(1+r)
    ## we might even not need growth in infected patches at all since we already have hosts movement to replenish infected patches
    ## just including this to ensure burnout can happen to patches that persisted after first wave

    ## update state
    N[t+1]=n
    P[t+1]=p
  }
  
  data.frame(time = seq.int(tt), N, P)
}


res <- simfun()
print(res)

## showing everything on one log plot works well
library(tinyplot)
res_long <- tidyr::pivot_longer(res, -time, names_to = "var")
par(las=1)
tinyplot(value ~ time | var, data = res_long, type = "l", log = "y")

## or more traditionally
plotfun2 <- function(res, ...) {
  par(mfrow=c(1,2))
  plot(P ~ time, data = res, type = "l", ...)
  matplot(res$time, res[c("N")], col = c(2,4), type = "l", ylab = "pop density", ...)
  legend("right", col = c(2,4), lty = 1:2, legend = c("N"))
}

plotfun2(res, lwd = 2)
