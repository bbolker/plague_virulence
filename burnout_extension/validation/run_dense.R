source('validation/R/scan.R')
Rvals<-c(1.05,seq(1.1,1.5,.1),1.75,2,2.5,3,3.5,4,5,6,7,10)
rhovals<-c(.0025,.005,.0075,.01,.02,.03,.05,.075,.1)
grid<-expand.grid(theta=c(0,.25,.5,.75,1),R0=Rvals,rho=rhovals,K=1e4)
run_grid(grid,'validation/data/dense_results.csv')

# K sensitivity subset
kg<-expand.grid(theta=c(0,.5,1),R0=c(1.2,2,3,5),rho=c(.005,.02,.05),K=c(1e3,1e4,1e5))
run_grid(kg,'validation/data/K_sensitivity.csv')
