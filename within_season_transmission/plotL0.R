library(burnout)

L0<-function(c0=0.5,
             S=1e5,
             R0=3,
             D=1,
             r=0.5,
             epsilon=0.02,
             a=0.3 #coloring everthing when setting a=0.5, but I don't know why
) {
  z=final_size(R0)
  N_DFE=c0*S/(1-(1+r)*(1-z*D)*(1-c0))
  P1<-P1_prob(R0,epsilon=epsilon,k=1,N=N_DFE)
  return(P1+a*z*(1-1/R0))
}


n <- 101
epsvec <- seq(1e-4, 0.05, length.out = n) ## warnings if eps=0
R0vec <- seq(1, 5, length.out = n)
grad_R0 <- L <- matrix(NA, length(epsvec), length(R0vec))
for (i in seq_along(epsvec)) {
  for (j in seq_along(R0vec)) {
    L[i,j] <- L0(epsilon = epsvec[i],R0 = R0vec[j]) 
    cat(" eps=", epsvec[i], " R0=", R0vec[j], " L0=", L[i,j], "\n")
  }
}
dR0 <- diff(R0vec)[1]
for (j in 2:length(R0vec)) {
  grad_R0[,j] <- (L[,j]-L[,j-1])/dR0
  cat(" eps=", epsvec[i], " R0=", R0vec[j], " grad_R0=", grad_R0[i,j], "\n")
}


par(las=1, bty = "l")
image(epsvec, R0vec, grad_R0>0, col = c(adjustcolor("blue", alpha = 0.3), "#FFFFFF"),
      main = "L0")
contour(epsvec, R0vec, L,
        levels = seq(0,2,by = 0.2),
        add = TRUE)
contour(epsvec, R0vec, grad_R0, level = 0, add = TRUE, col = 2, lwd = 2)