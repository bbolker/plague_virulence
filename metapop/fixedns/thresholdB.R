library(burnout)

getDFE<-function(c0=0.5,
                 S=1e6,
                 R0=3,
                 D=1,
                 r=0.5
) {
  
  z=final_size(R0)
  N_DFE=c0*S/(1-(1+r)*(1-z*D)*(1-c0))
  return(N_DFE)}


P1_at_DFE<-function(R0=2,epsilon=0.01){
  N_DFE<-getDFE(R0=R0)
  P1<-P1_prob(R0,epsilon=epsilon,k=1,N=N_DFE)
}

thresholdB<-function(R0=2,epsilon=0.01){
  N_DFE<-getDFE(R0=R0)
  P1<-P1_prob(R0,epsilon=epsilon,k=1,N=N_DFE)
  B=1/P1-1
}

n <- 101
epsvec <- seq(1e-4, 0.02, length.out = n) ## warnings if eps=0
R0vec <- seq(1, 2.5, length.out = n)
B<-grad_R0 <- persist <- matrix(NA, length(epsvec), length(R0vec))
for (i in seq_along(epsvec)) {
  for (j in seq_along(R0vec)) {
    persist[i,j] <- P1_at_DFE(R0 = R0vec[j], epsilon = epsvec[i]) 
    B[i,j]<- thresholdB(R0 = R0vec[j], epsilon = epsvec[i])
  }
}
dR0 <- diff(R0vec)[1]
for (j in 2:length(R0vec)) {
  grad_R0[,j] <- (persist[,j]-persist[,j-1])/dR0
}


par(las=1, bty = "l")

image(epsvec, R0vec, grad_R0>0, col = c(adjustcolor("blue", alpha = 0.3), "#FFFFFF"),
      main = "threshold logB")

contour(epsvec, R0vec, log10(B),
        levels = c(-1,0,0.5,1,2,4,8),
        add = TRUE)
contour(epsvec, R0vec, grad_R0, level = 0, add = TRUE, col = 2, lwd = 2)

