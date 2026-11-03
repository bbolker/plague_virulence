source('validation/R/stochastic.R');source('validation/R/native_simulation.R')
outdir<-'validation/data/tau_bias_audit_20261104'
exact<-compile_native_simulation('validation/src/exact_first_trough.c')
guard<-compile_native_simulation('validation/src/guarded_first_trough.c')
cl<-parallel::makeCluster(4L)
parallel::clusterCall(cl,function(p){dyn.load(p[1]);dyn.load(p[2]);NULL},c(exact,guard))
design<-data.frame(method='native_guard',event_budget=c(.1,.03,.01,.003),
 epsilon=.01,small_count=100,tf=500,R0=128,rho=.01,theta=.25,K=1e6,n=10000L,
 seed=202611800L+1:4)
write.csv(design,file.path(outdir,'native_guard_design.csv'),row.names=FALSE)
results<-list()
for(i in seq_len(nrow(design))){
 row<-design[i,];tot<-NULL;start<-proc.time()[3]
 for(batch in 1:20){
  z<-native_batch(cl,row,500,row$seed+batch*100L)
  tot<-if(is.null(tot))z else tot+z
  cat('budget',row$event_budget,'attempts',tot['attempts'],'\n');flush.console()
 }
 ci<-wilson(tot['persistent'],tot['established'])
 results[[i]]<-cbind(row,as.data.frame(as.list(tot)),data.frame(P=tot['persistent']/tot['established'],
   CI_low=max(0,ci[1]),CI_high=min(1,ci[2]),seconds=proc.time()[3]-start))
 write.csv(do.call(rbind,results),file.path(outdir,'native_guard_results.csv'),row.names=FALSE)
 print(results[[i]],row.names=FALSE);stopifnot(tot['unresolved']==0)
}
parallel::stopCluster(cl)
