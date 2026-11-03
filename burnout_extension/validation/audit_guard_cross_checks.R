source('validation/R/stochastic.R');source('validation/R/native_simulation.R')
outdir<-'validation/data/tau_bias_audit_20261104'
exact<-compile_native_simulation('validation/src/exact_first_trough.c')
guard<-compile_native_simulation('validation/src/guarded_first_trough.c')
cl<-parallel::makeCluster(8L)
parallel::clusterCall(cl,function(p){dyn.load(p[1]);dyn.load(p[2]);NULL},c(exact,guard))
# Parameter choices and sample sizes fixed before inspecting this audit.
base<-data.frame(R0=c(16,32,64,96,1.05,1.05),theta=c(.25,.25,.25,.25,0,.5),
 K=c(1e6,1e5,1e6,1e6,1e6,1e6))
design<-base[rep(seq_len(nrow(base)),each=2),]
design$method<-rep(c('exact','native_guard'),nrow(base))
design$rho<-.01;design$epsilon<-.01;design$event_budget<-.01
design$small_count<-100;design$tf<-2000;design$n<-4000L
# Additional independent large-K step-refinement checks.
big<-data.frame(R0=rep(c(96,128,96,128),each=2),theta=.5,
 K=rep(c(1e8,1e8,1e9,1e9),each=2),method='native_guard',rho=.01,
 epsilon=rep(c(.01,.003),4),event_budget=rep(c(.01,.003),4),
 small_count=100,tf=500,n=4000L)
design<-rbind(design,big);design$seed<-202611900L+seq_len(nrow(design))
write.csv(design,file.path(outdir,'cross_check_design.csv'),row.names=FALSE)
results<-list()
for(i in seq_len(nrow(design))){
 row<-design[i,];tot<-NULL;start<-proc.time()[3]
 for(batch in 1:8){
  z<-native_batch(cl,row,500,row$seed+batch*100L)
  tot<-if(is.null(tot))z else tot+z
 }
 ci<-wilson(tot['persistent'],tot['established'])
 results[[i]]<-cbind(row,as.data.frame(as.list(tot)),data.frame(P=tot['persistent']/tot['established'],
   CI_low=max(0,ci[1]),CI_high=min(1,ci[2]),seconds=proc.time()[3]-start))
 write.csv(do.call(rbind,results),file.path(outdir,'cross_check_results.csv'),row.names=FALSE)
 print(results[[i]],row.names=FALSE);flush.console();stopifnot(tot['unresolved']==0)
}
parallel::stopCluster(cl)
