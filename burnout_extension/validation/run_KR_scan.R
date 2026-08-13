source('validation/R/scan.R')
# Uniform in log(R0-1), as used by the requested horizontal plot coordinate.
R0vals<-1+exp(seq(log(.05),log(9),length.out=15))
grid<-expand.grid(theta=c(0,.25,.5,.75,1),R0=R0vals,
                  rho=c(.02,.04,.06,.08,.10),
                  K=c(1e3,3e3,1e4,3e4,1e5))
run_grid(grid,'validation/data/K_R0_scan.csv')
