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
#' @param N number of patches
simfun <- function(tt  = 100,
                   B0 = 0.8,
                   r = 0.5,
                   D = 0.5,
                   epsilon=0.05,
                   R0=10,
                   ki=0.01,
                   ks=0.01,
                   c0=0.2,
                   K=300,
                   N = 1000) {

  if (!require("burnout")) stop("please install the 'burnout' package: ",
                                "`remotes::install_github('davidearn/burnout')`")

  ## Initialize parameters
  P <- NI <- NS <- rep(NA, tt)

  ## Initial conditions
  ## FIXME: don't hardcode?
  P[1] = 0.01                  
  NI[1] = 100
  NS[1] = 100

### Simulation loop
  for (t in 1:(tt-1)) {
    ## Current state
    p = P[t]
    ni = NI[t]
    ns = NS[t]
    
    if (!is.na(p) && p<1/N) {
      ## 'stochastic' (i.e. discrete) extinction
      p=0
      ni=0
      break
    } else{
      ## migration and infection
      B=B0*(1-exp(-ki*ni))*(1-exp(-ks*ns))
      delta_p=B*p*(1-p)
      ni=(p*ni+delta_p*ns)/(p+delta_p)
      p = p + delta_p 

      c=min(c0,1/(1/p+1/(1-p)))
      delta_N <- c*(ns-ni)
      ni <- ni + delta_N/p
      ns <- ns - delta_N/(1-p)
      
      ##fizzle and burnout
      P1=P1_prob(R0,epsilon,k=1,N=ni)
      delta_p<-(1-P1)*p
      ns<-((1-p)*ns+delta_p*ni)/(1-p+delta_p)
      p <-p-delta_p
      
      ##epidemic
      z=final_size(R0)
      ni <- ni * (1 - z*D)  
    }
    
    ## birth
    ns=ns+r*ns*(1-ns/K)
    ni=ni+r*ni*(1-ni/K)
    
    NI[t+1]=ni
    NS[t+1]=ns
    P[t+1]=p
  }

  data.frame(time = seq.int(tt), NI, NS, P)
}


res <- simfun()

## showing everything on one log plot works well
library(tinyplot)
res_long <- tidyr::pivot_longer(res, -time, names_to = "var")
par(las=1)
tinyplot(value ~ time | var, data = res_long, type = "l", log = "y")

## or more traditionally
par(mfrow=c(1,2))
plot(P ~ time, data = res, type = "l")
matplot(res$time, res[c("NI", "NS")], col = c(2,4), type = "l", ylab = "pop density")
legend("right", col = c(2,4), lty = 1:2, legend = c("NI", "NS"))
