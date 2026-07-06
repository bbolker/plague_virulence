library("burnout")

# Function to compute DFE
getDFE <- function(c0=0.5,
                   S=1e6,
                   R0=3,
                   D=1,
                   r=0.5) {
  z <- final_size(R0)
  N_DFE <- c0*S/(1-(1+r)*(1-z*D)*(1-c0))
  return(N_DFE)
}

# Function to compute P1 at DFE
P1_at_DFE <- function(R0=2, epsilon=0.02, c0=0.5, r=0.5){
  N_DFE <- getDFE(R0=R0, c0=c0, r=r)
  P1 <- P1_prob(R0, epsilon=epsilon, k=1, N=N_DFE)
  return(P1)
}

# Parameter grid
n <- 101
epsvec <- seq(1e-4, 0.03, length.out = n)
R0vec <- seq(1, 2.5, length.out = n)

# Save all plots to one pdf
pdf("more_P1atDFE.pdf", width=8, height=8)

for (c0_val in c(0.25, 0.5, 1)) {
  for (r_val in c(0.25, 0.5, 1)) {
    
    # Compute persistence matrix
    persist <- matrix(NA, length(epsvec), length(R0vec))
    grad_R0 <- matrix(NA, length(epsvec), length(R0vec))
    
    for (i in seq_along(epsvec)) {
      for (j in seq_along(R0vec)) {
        persist[i,j] <- P1_at_DFE(R0=R0vec[j], epsilon=epsvec[i],
                                  c0=c0_val, r=r_val)
      }
    }
    
    dR0 <- diff(R0vec)[1]
    for (j in 2:length(R0vec)) {
      grad_R0[,j] <- (persist[,j] - persist[,j-1]) / dR0
    }
    
    # Plot figure
    par(las=1, bty="l", mar=c(4,4,3,2))
    image(epsvec, R0vec, grad_R0 > 0,
          col=c(adjustcolor("blue", alpha=0.3), "#FFFFFF"),
          xlab="epsilon", ylab="R0",
          main=paste0("P1 at DFE (c0=",c0_val,", r=",r_val,")"))
    
    contour(epsvec, R0vec, log10(persist),
            levels=c(-12, -8, -4, -2, -1, log10(0.2), log10(0.3)),
            add=TRUE)
    contour(epsvec, R0vec, grad_R0, level=0, add=TRUE, col=2, lwd=2)
  }
}

dev.off()
