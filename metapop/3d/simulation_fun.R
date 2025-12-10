## avoid T, it's a built-in synonym for TRUE
#' @param tt number of time steps
#' @param alpha infection rate between patches
#' @param r  growth rate of rat population
#' @param D  death rate due to pathogen
#' @param epsilon ratio of pathogen to host generation time
#' @param R0 intrinsic reproductive number of pathogen
#' @param c0
#' @param K  rat carrying capacity per patch
#' @param N number of patches (only affects discrete-extinction logic)
#' @param starting conditions

simfun <- function(tt  = 300,
                   alpha=5e-6,
                   r = 0.5,
                   D = 1,
                   epsilon=0.05,
                   R0=2.5,
                   c0=0.01,
                   K=1e6,
                   n = 100,
                   start = c(p=0.01, NS=K, NI=K)) {
  
  if (!require("burnout")) stop("please install the 'burnout' package: ",
                                "`remotes::install_github('davidearn/burnout')`")
  
  ## Initialize parameters
  P <- NI <- NS <- rep(NA, tt)
  
  P[1] <- start[["p"]]
  NS[1] <- start[["NS"]]
  NI[1] <- start[["NI"]]
  
  ## Simulation loop
  for (t in 1:(tt-1)) {
    ## Current state
    p = P[t]
    ni = NI[t]
    ns = NS[t]

    # if (!is.na(p) && p<1/N) {
    #   ## 'stochastic' (i.e. discrete) extinction
    #   p=0
    #   ni=0
    #   break
    # } else{
     if (p==0) {
      ## continuous version, derived from the limit when p approaches 0

      ## burnout
      ni <- ni * (1 - final_size(R0)*D)
      
      ## logistic growth
      ni <- ni + r*ni*(1-ni/K)
      ns <- ns + r*ns*(1-ns/K)
       
      ## host movement
      ni <- ni + c0*(ns-ni)

      ## infection between patches
      B=alpha*c0*n*p*z*NI[t]
      ni<-(ni+B*ns)/(1+B)
      
    } else{
      
      ## burnout
      Pb <- burnout_prob(R0,epsilon,N=ni)
      z=final_size(R0)
      ni=ni*(1-z*D)
      ns=(ni*Pb*p+ns*(1-p))/(Pb*p+1-p)
      p <- (1-Pb)*p
      
      ## birth
      ns=ns+r*ns*(1-ns/K)
      ni=ni+r*ni*(1-ni/K)
      
      ## host movement between patches
      delta_N <- c0*(ns-ni)
      ni <- ni + delta_N*(1-p)
      ns <- ns - delta_N*p
      
      ## infection between patches
      delta_p=(1-exp(-alpha*c0*n*p*z*NI[t]))*(1-p)    ## infection of S patches
      ni=(p*ni+delta_p*ns)/(p+delta_p)  ## ni changes because S patches convert to I
      p = p + delta_p 
      
      
      ##fizzle
      Pf=1/R0
      ns=((1-p)*ns+Pf*p*ni)/(1-p+Pf*p)
      p=(1-Pf)*p
      
    }
    
    NI[t+1]=ni
    NS[t+1]=ns
    P[t+1]=p
  }
  
  data.frame(time = seq.int(tt), NI, NS, P)
}

res <- simfun(R0=2.5)

## showing everything on one log plot works well
library(tinyplot)
res_long <- tidyr::pivot_longer(res, -time, names_to = "var")
par(las=1)
tinyplot(value ~ time | var, data = res_long, type = "l", log = "y")

## or more traditionally
plotfun2 <- function(res, ...) {
  par(mfrow=c(1,2))
  plot(P ~ time, data = res, type = "l", ...)
  matplot(res$time, res[c("NI", "NS")], col = c(2,4), type = "l", ylab = "pop density", ...)
  legend("right", col = c(2,4), lty = 1:2, legend = c("NI", "NS"))
}

#plotfun2(res, lwd = 2)


