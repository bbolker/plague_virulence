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

simfun <- function(tt  = 200,
                   c0=1,
                   S = 300,
                   B0 = 1.5,
                   R0 = 2.5,
                   epsilon=0.05,
                   #D=0.85,
                   D=1,
                   r = 0.5,
                   K = 300,
                   start = c(p=0.01,N=S)) {
  
  if (!require("burnout")) stop(
    "please install the 'burnout' package: ",
    "`remotes::install_github('davidearn/burnout')`")
  
  z=final_size(R0)
  if((1+r)*(1-z*D)>1){
    cat("growth rate too large/death rate too low\n")
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
    
    ##epidemic:death of hosts and burnout 
    ##burnout
    P1=P1_prob(R0,epsilon,k=1,N=n)
    n <- n * (1 - p*z*D*(1-1/R0))*(1+r-r*n/K)  ##death of hosts
    p <-P1*p+(1-exp(-B0*p))*(1-p)  ##infection (hazard model)
    
    ## update state
    N[t+1]=n
    P[t+1]=p
  }
  
  data.frame(time = seq.int(tt), N, P)
}


res1 <- simfun(D=1)
#res2 <- simfun(D=0.85)

## showing everything on one log plot works well
library(tinyplot)
res_long <- tidyr::pivot_longer(res1, -time, names_to = "var")
par(las=1)
tinyplot(value ~ time | var, data = res_long, type = "l", log = "y")

## or more traditionally
plotfun2 <- function(res1,...) {
  par(mfrow=c(1,2))
  plot(res1$time, res1$P,
          type = "l", col = c("red", "blue"), lty = 1, lwd = 2,
          xlab = "Time", ylab = "Fraction infected")
  
  
  plot(res1$time, res1$N,
          type = "l", col = c("red", "blue"), lty = 1, lwd = 2,
          xlab = "Time", ylab = "pop density")

}

plotfun2(res1,lwd = 2)

## plot trajectories
plot_traj <- function(
    p_seq = seq(0, 1, by = 0.2),
    N_seq = seq(0, 300, by = 100),
    ...) {
  par(mfrow=c(1,1))
  all_traj <- list()
  traj_labels <- c()
  
  for (p0 in p_seq) {
    for (N0 in N_seq) {
      traj <- simfun(start = c(p = p0, N = N0))
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
