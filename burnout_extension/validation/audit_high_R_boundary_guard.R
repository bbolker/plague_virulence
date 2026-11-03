# Small fixed-design audit before any production-data replacement.
source('validation/R/stochastic.R');source('validation/R/tau_boundary_guard.R')
outdir<-'validation/data/tau_bias_audit_20261104'
cl<-parallel::makeCluster(4L)
parallel::clusterExport(cl,c('one_boundary_guarded_tau','tau_boundary_max_step','classify_first_trough'))
design<-data.frame(event_budget=c(.1,.03,.01),epsilon=.01,n=4000,
                  seed=202611600L+1:3,R0=128,rho=.01,theta=.25,K=1e6,tf=100)
write.csv(design,file.path(outdir,'guard_design.csv'),row.names=FALSE)
results<-list()
for(i in seq_len(nrow(design))){
 row<-design[i,];totals<-c(attempts=0,established=0,persistent=0,unresolved=0,ever_zero=0)
 start<-proc.time()[3]
 for(batch in 1:8){
  parallel::clusterSetRNGStream(cl,iseed=row$seed+batch*100L)
  a<-parallel::parLapply(cl,rep(125L,4),function(n,row){
   z<-rowSums(replicate(n,one_boundary_guarded_tau(row$R0,row$rho,row$theta,row$K,
                row$epsilon,row$event_budget,tf=row$tf)))
   c(attempts=n,established=unname(z[1]),persistent=unname(z[2]),
     unresolved=unname(z[3]),ever_zero=unname(z[4]))
  },row)
  totals<-totals+colSums(do.call(rbind,a))
  cat('budget',row$event_budget,'attempts',totals['attempts'],'\n');flush.console()
 }
 ci<-wilson(totals['persistent'],totals['established'])
 results[[i]]<-cbind(row,data.frame(established=totals['established'],persistent=totals['persistent'],
  unresolved=totals['unresolved'],ever_zero=totals['ever_zero'],
  P=totals['persistent']/totals['established'],CI_low=max(0,ci[1]),CI_high=min(1,ci[2]),
  seconds=proc.time()[3]-start))
 write.csv(do.call(rbind,results),file.path(outdir,'guard_results.csv'),row.names=FALSE)
 print(results[[i]],row.names=FALSE);stopifnot(totals['unresolved']==0)
}
parallel::stopCluster(cl)
