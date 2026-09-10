source('2d_ode_validation/R/trough_kendall_psi.R')
stopifnot(requireNamespace('data.table',quietly=TRUE))
library(data.table)

sim<-fread('psi_validation/data/psi05_stochastic_results.csv')
dense<-CJ(rho=c(.01,.02,.05,.10),theta=c(0,.5,1),
  K=c(1000,3000,10000,30000),
  R0=1+exp(seq(log(.05),log(5),length.out=101)))
dense[,psi:=.5]
dense[,R0:=round(R0,12)]
sim[,R0:=round(R0,12)]
grid<-unique(rbindlist(list(dense,sim[,.(rho,theta,K,R0,psi)])))
setorder(grid,rho,theta,K,R0)
ans<-vector('list',nrow(grid))
for(i in seq_len(nrow(grid))){
  g<-grid[i]
  ans[[i]]<-tryCatch(
    trough_kendall_row_2d(g$R0,g$rho,g$theta,g$psi,g$K),
    error=function(e)data.frame(rho=g$rho,theta=g$theta,K=g$K,R0=g$R0,
      psi=g$psi,status=paste0('ERROR: ',conditionMessage(e))))
  if(i%%50L==0L)cat(sprintf('%d/%d deterministic trajectories complete\n',i,nrow(grid)))
}
curves<-rbindlist(ans,fill=TRUE)
curves[,`:=`(P_trough_unconditional=(1-1/R0)*P_trough,
  P_kendall_unconditional=(1-1/R0)*P_kendall)]
dir.create('2d_ode_validation/data',FALSE,TRUE)
fwrite(curves,'2d_ode_validation/data/psi05_2d_ode_curves.csv')
cat('2D ODE and Kendall curves complete\n')
