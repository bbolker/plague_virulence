getDFE<-function(c0=0.5,
                 S=1e6,
                 R0=3,
                 D=1,
                 r=0.5
) {
  if (!require("burnout")) stop(
    "please install the 'burnout' package: ",
    "`remotes::install_github('davidearn/burnout')`")
  
  z=final_size(R0)
  N_DFE=c0*S/(1-(1+r)*(1-z*D)*(1-c0))
  return(N_DFE)}


P1_at_DFE<-function(R0=2,epsilon=0.01){
  N_DFE<-getDFE(R0=R0)
  P1<-P1_prob(R0,epsilon=epsilon,k=1,N=N_DFE)
}


plot_P1 <- function(P_func=P1_at_DFE, x_start = 1, x_end = 5, step = 0.01) {
  x_vals <- seq(x_start, x_end, by = step)
  y_vals <- P_func(x_vals)
  
  plot(x_vals, y_vals, type = "l", col = "blue", lwd = 2,
       xlab = "R0", ylab = "Persistence prob", main = "change of P1 with R0")
}

plot_P1()
