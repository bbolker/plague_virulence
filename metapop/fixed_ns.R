## avoid T, it's a built-in synonym for TRUE
#' @param tt number of time steps
#' @param B0 infection rate between patches
#' @param r  growth rate of rat population
#' @param D  death rate due to pathogen
#' @param epsilon ratio of pathogen to host generation time
#' @param R0 intrinsic reproductive number of pathogen
#' @param c0
#' @param K  rat carrying capacity per patch
#' @param starting conditions


getDFE<-function(c0=0.2,
                 S=300,
                 R0=10,
                 D=0.5,
                 r=0.5
) {
  if (!require("burnout")) stop(
    "please install the 'burnout' package: ",
    "`remotes::install_github('davidearn/burnout')`")
  
  z=final_size(R0)
  N_DFE=c0*(1+r)*(1-z*D)*S/(1-(1+r)*(1-z*D)*(1-c0))
  return(N_DFE)}

canPersist<-function(c0=0.2,
                     S=300,
                     B0=0.5,
                     R0=10,
                     epsilon=0.05,
                     D=0.5,
                     r=0.5){
  N_DFE<-getDFE(c0=c0,S=S,R0=R0,D=D,r=r)
  P1<-P1_prob(R0,epsilon,k=1,N=N_DFE+c0*(S-N_DFE))
  eigenvalue<-(1+B0)*P1
  return (eigenvalue>=1) 
}

getEE<-function(c0=0.2,
                S = 300,
                B0 = 0.8,
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
  
  z=final_size(R0)
  if((1+r)*(1-z*D)>1){
    cat("growth rate too large/death rate too low")
    return(invisible(NULL))}
  
  if(!canPersist(c0=c0,
                 S=S,
                 B0=B0,
                 R0=R0,
                 epsilon=epsilon,
                 D=D,
                 r=r)){
    cat("epidemic cannot persist")
    return(invisible(NULL))}
  
  f <- function(x) {
    p<-x[1]
    n<-x[2]
    
    P1=P1_prob(R0,epsilon,k=1,N=n+c0*(1-p)*(S-n))
    
    eq1<-P1*(1+B0*(1-p))-1
    eq2<-(1+r)*(1-z*D)*(n+c0*(1-p)*(S-n))-n
    
    return(c(eq1, eq2))
  }
  
  x0 <- c(0,8)  
  result <- nleqslv(x0, f)
  
  return(result$x)
}

simfun <- function(tt  = 100,
                   c0=0.2,
                   S = 300,
                   B0 = 0.8,
                   R0=10,
                   epsilon=0.05,
                   D = 0.5,
                   r = 0.5,
                   start = c(p=0.01,N=S)) {
  
  if (!require("burnout")) stop(
    "please install the 'burnout' package: ",
    "`remotes::install_github('davidearn/burnout')`")
  
  
  z=final_size(R0)
  if((1+r)*(1-z*D)>1){
    cat("growth rate too large/death rate too low")
    return(invisible(NULL))}
  
  ## Initialize parameters
  P <- N <- rep(NA, tt)  ## P for fraction of infected patches, N for population density of infected patches
  
  P[1] <- start[["p"]]
  N[1] <- start[["N"]]
  
  ## Simulation loop
  for (t in 1:(tt-1)) {
    ## Current state
    p <- P[t]
    n <- N[t]
    
    ##ignore stocastic extinction
    
    ## host movement and infection between patches(caused by host movement)
    n <- n + c0*(S-n)*(1-p)            ## host (rat) movement between patches
    p <- p + B0*p*(1-p)                 ##infection of S patches (constant B for simplicity)
    ## ignore ni change of n because S patches convert to I (combined with movement event for simplicity)
    
    ##epidemic:death of hosts and burnout 
    ##burnout
    P1=P1_prob(R0,epsilon,k=1,N=n)
    p <-P1*p 
    ##death of hosts
    n <- n * (1 - z*D) 
    
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

## plot trajectories
plot_traj <- function(B0=0.8,R0=10,r=0.5,
                      p_seq = seq(0, 1, by = 0.2),
                      N_seq = seq(0, 300, by = 100),
                      steps = 100,
                      ...) {
  par(mfrow=c(1,1))
  all_traj <- list()
  traj_labels <- c()
  
  for (p0 in p_seq) {
    for (N0 in N_seq) {
      traj <- simfun(B0 = B0, R0 = R0, r = r, start = c(p = p0, N = N0), tt = steps)
      all_traj[[length(all_traj) + 1]] <- traj
    }
  }
  
  P_list <- lapply(all_traj, function(df) df$P)
  N_list <- lapply(all_traj, function(df) df$N)
  
  max_len <- max(sapply(P_list, length))
  P_mat <- do.call(cbind, lapply(P_list, function(x) c(x, rep(NA, max_len - length(x)))))
  N_mat <- do.call(cbind, lapply(N_list, function(x) c(x, rep(NA, max_len - length(x)))))
  
  matplot(P_mat, N_mat, type = "l", lty = 1, lwd=2, col = rainbow(length(all_traj)),
          xlab = "p",
          ylab = "N",
          main = "Trajectories", ...)
  }
plot_traj()
