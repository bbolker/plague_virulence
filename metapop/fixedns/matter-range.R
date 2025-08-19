library("burnout")

canPersist<-function(c0=0.5,
                     S=300,
                     B0=10,
                     R0=2.5,
                     epsilon=0.05,
                     D=1,
                     r=0.5){
  z=final_size(R0)
  N_DFE=c0*S/(1-(1+r)*(1-z*D)*(1-c0))
  P1<-P1_prob(R0,epsilon,k=1,N=N_DFE)
  eigenvalue<-(1+B0)*P1
  return (eigenvalue>=1) 
}

matter<-function(epsilon=0.05,
                 R0=3){
    return(canPersist(D=1,R0=R0,epsilon=epsilon)!=canPersist(D=0.85,R0=R0,epsilon=epsilon))
  }

n <- 101
epsvec <- seq(1e-4, 0.2, length.out = n) ## warnings if eps=0
R0vec <- seq(1, 3, length.out = n)
M <- matrix(NA, length(epsvec), length(R0vec))
for (i in seq_along(epsvec)) {
  for (j in seq_along(R0vec)) {
    M[i,j] <- matter(R0=R0vec[j], epsilon=epsvec[i])
  }
}

image(epsvec, R0vec, M, col=c("#FFFFFF", adjustcolor("blue", alpha=0.3)))



