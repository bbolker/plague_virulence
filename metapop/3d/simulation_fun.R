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

getDFE<-function(c0=0.2,
                 K = 300,
                 B0 = 0.5,
                 R0=10,
                 epsilon=0.05,
                 D = 1,
                 r = 0.5
) {
  if (!require("burnout")) stop(
    "please install the 'burnout' package: ",
    "`remotes::install_github('davidearn/burnout')`")
  
  if (!require("nleqslv")) stop(
    "Please install the 'nleqslv' package:\n",
    "install.packages('nleqslv')")
  
  f <- function(x) {
    ni<-x
    ni1 <- ni + c0*(K-ni)
    B=B0
    ni2<-(ni1+B*K)/(1+B)
    ni3 <- ni2 * (1 - final_size(R0)*D)
    ni4 <- ni3 + r*ni3*(1-ni3/K)
    eq1<-ni4-ni
    return(c(eq1))
  }
  
  x0 <- c(1)  
  result <- nleqslv(x0, f)
  
  return(result$x)
}


getEE<-function(c0=0.2,
                K = 300,
                B0 = 0.5,
                R0=10,
                epsilon=0.05,
                D = 0.5,
                r = 0.5
) {
  if (!require("burnout")) stop(
    "please install the 'burnout' package: ",
    "`remotes::install_github('davidearn/burnout')`")
  
  if (!require("nleqslv")) stop(
    "Please install the 'nleqslv' package:\n",
    "install.packages('nleqslv')")
  
  f <- function(x) {
    p<-x[1]
    ni<-x[2] 
    ns<-x[3]
    
    delta_N <- c0*(ns-ni)
    p1<-p
    ni1 <- ni + delta_N*(1-p)
    ns1 <- ns - delta_N*p
    
    B=B0
    p2 = p1 + B*p1*(1-p1)  
    ni2=(ni1+B*(1-p1)*ns1)/(1+B*(1-p1))
    ns2<-ns1
    
    P1=P1_prob(R0,epsilon,k=1,N=ni2)
    delta_p<-(1-P1)*p2
    p3 <-p2-delta_p
    ni3<-ni2
    ns3<-((1-p2)*ns2+delta_p*ni2)/(1-p2+delta_p)
    
    z=final_size(R0)
    p4<-p3
    ni4 <- ni3 * (1 - z*D)  
    ns4<-ns3
    
    ns5=ns4+r*ns4*(1-ns4/K)
    ni5=ni4+r*ni4*(1-ni4/K)
    p5<-p4
    
    eq1<-p5-p
    eq2<-ni5-ni
    eq3<-ns5-ns
    
    return(c(eq1, eq2,eq3))
  }
  
  x0 <- c(0.2,100,200)  
  result <- nleqslv(x0, f)
  
  return(result$x)
}


simfun <- function(tt  = 300,
                   B0 = 0.5,
                   r = 0.5,
                   D = 1,
                   epsilon=0.05,
                   R0=2.5,
                   ki=0.01,
                   ks=0.01,
                   c0=0.2,
                   K=300,
                   N = 1000,
                   start = c(p=0, NS=K, NI=K)) {
  
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
      
      ## host movement
      ni <- ni + c0*(ns-ni)
      
      ## infection between patches
      B=B0
      ni<-(ni+B*ns)/(1+B)
      
      ## epidemic
      ni <- ni * (1 - final_size(R0)*D)
      
      ## logistic growth
      ni <- ni + r*ni*(1-ni/K)
      ns <- ns + r*ns*(1-ns/K)
      
    } else{
      ## migration and infection
      
      ## host (rat) movement between patches
      delta_N <- c0*(ns-ni)
      ni <- ni + delta_N*(1-p)
      ns <- ns - delta_N*p
      
      ## infection between patches
      #B=B0*(1-exp(-ki*ni))*(1-exp(-ks*ns))
      B=B0
      delta_p=B*p*(1-p)                 ## infection of S patches
      ni=(ni+B*(1-p)*ns)/(1+B*(1-p))  ## ni changes because S patches convert to I
      p = p + delta_p 
      
      
      ##fizzle and burnout
      P1=P1_prob(R0,epsilon,k=1,N=ni)
      delta_p<-(1-P1)*p
      ns<-((1-p)*ns+delta_p*ni)/(1-p+delta_p)
      p <-p-delta_p
      
      ##epidemic
      z=final_size(R0)
      ni <- ni * (1 - z*D)  
      
      ## birth
      ns=ns+r*ns*(1-ns/K)
      ni=ni+r*ni*(1-ni/K)
    }
    
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
plotfun2 <- function(res, ...) {
  par(mfrow=c(1,2))
  plot(P ~ time, data = res, type = "l", ...)
  matplot(res$time, res[c("NI", "NS")], col = c(2,4), type = "l", ylab = "pop density", ...)
  legend("right", col = c(2,4), lty = 1:2, legend = c("NI", "NS"))
}

plotfun2(res, lwd = 2)

res3 <-  simfun(start = c(p=0.1, NS=300, NI=300))
plotfun2(res3, lwd = 2)




