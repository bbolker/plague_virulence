library("burnout")

# Parameter sets
S_vals <- round(exp(seq(log(1e2), log(1e6), length.out = 5)))  # 5 values from 1e2 to 1e6
B0_vals <- c(1,5,10,15,20)                          # 5 values from 1 to 20

# Set 5x5 plotting layout
par(mfrow=c(5,5), mar=c(3,3,2,1), oma=c(2,2,2,2))  # mar: inner margins, oma: outer margins

n <- 101
epsvec <- seq(1e-4, 0.2, length.out = n)
R0vec <- seq(1, 3, length.out = n)

for (s in S_vals) {
  for (b in B0_vals) {
    
    # Current canPersist function with specific S and B0
    canPersist_curr <- function(c0=0.5,
                                S=s,
                                B0=b,
                                R0=2.5,
                                epsilon=0.05,
                                D=1,
                                r=0.5){
      z <- final_size(R0)
      N_DFE <- c0*S/(1-(1+r)*(1-z*D)*(1-c0))
      P1 <- P1_prob(R0, epsilon, k=1, N=N_DFE)
      eigenvalue <- (1+B0)*P1
      return(eigenvalue >= 1) 
    }
    
    # Current matter function
    matter_curr <- function(epsilon=0.05, R0=3){
      return(canPersist_curr(D=1,R0=R0,epsilon=epsilon) != canPersist_curr(D=0.85,R0=R0,epsilon=epsilon))
    }
    
    # Compute matrix
    M <- matrix(NA, length(epsvec), length(R0vec))
    for (i in seq_along(epsvec)) {
      for (j in seq_along(R0vec)) {
        M[i,j] <- matter_curr(R0=R0vec[j], epsilon=epsvec[i])
      }
    }
    
    # Plot the matrix
    image(epsvec, R0vec, M, col=c("#FFFFFF", adjustcolor("blue", alpha=0.3)),
          xlab="epsilon", ylab="R0", main=paste0("S=",s,", B0=",b))
  }
}
