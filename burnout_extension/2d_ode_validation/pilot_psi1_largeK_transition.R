source('psi_validation/R/stochastic_psi.R')
stopifnot(requireNamespace('data.table',quietly=TRUE))
library(data.table)

grid<-CJ(rho=c(.01,.02,.05,.10),K=c(1e6,1e7,1e8,1e9),
  q=c(2.5,3,3.5,4,4.5))
grid[,`:=`(theta=.5,psi=1,R0=rho*sqrt(K)/q,
  trough_scale=q*rho*sqrt(K))]
grid<-grid[R0>1.05]
setorder(grid,rho,K,R0)

cl<-parallel::makeCluster(min(16L,parallel::detectCores()))
on.exit(parallel::stopCluster(cl),add=TRUE)
parallel::clusterEvalQ(cl,library(adaptivetau))
parallel::clusterExport(cl,'one_adaptive_tau_psi_chunked',envir=environment())

outfile<-'2d_ode_validation/data/psi1_largeK_transition_pilot.csv'
out<-if(file.exists(outfile))fread(outfile)else data.table()
for(i in seq_len(nrow(grid))){
  g<-grid[i]
  if(nrow(out)&&out[rho==g$rho&K==g$K&q==g$q,.N])next
  n_pilot<-c('1e+06'=200L,'1e+07'=150L,'1e+08'=100L,'1e+09'=50L)[format(g$K,scientific=TRUE)]
  z<-simulate_counts_batch_tau_psi_chunked_lb(cl,g$R0,g$rho,g$theta,g$psi,g$K,
    n_attempts=n_pilot,seed=920000L+i*1000L,tau_epsilon=.01,chunk_time=50)
  row<-data.table(rho=g$rho,K=g$K,q=g$q,R0=g$R0,
    trough_scale=g$trough_scale,theta=g$theta,psi=g$psi,
    attempts=unname(z['attempts']),established=unname(z['established']),
    persistent=unname(z['persistent']),unresolved=unname(z['unresolved']))
  row[,P_conditional:=persistent/established]
  out<-rbind(out,row,fill=TRUE);setorder(out,rho,K,q);fwrite(out,outfile)
  cat(sprintf('%d/%d rho=%g K=%.0e R0=%.5g B=%g cond=%.3f unresolved=%d\n',
    i,nrow(grid),g$rho,g$K,g$R0,g$trough_scale,
    if(z['established']>0)z['persistent']/z['established'] else NA_real_,z['unresolved']))
}
