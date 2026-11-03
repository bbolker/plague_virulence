source('validation/R/native_simulation.R')
dll<-compile_native_simulation('validation/src/guarded_first_trough.c')
cl<-parallel::makeCluster(2L);parallel::clusterCall(cl,function(p)dyn.load(p),dll)
grid<-expand.grid(R0=c(8,32,128),theta=c(0,.25,.5,.75,1),K=c(1e6,1e7,1e8,1e9))
grid$method<-'native_guard';grid$rho<-.01;grid$epsilon<-.01
grid$event_budget<-.01;grid$small_count<-100;grid$tf<-500
out<-list()
for(i in seq_len(nrow(grid))){
 t<-proc.time()[3];z<-native_batch(cl,grid[i,],32,202613000L+i)
 out[[i]]<-cbind(grid[i,],data.frame(n=32,seconds=proc.time()[3]-t,unresolved=z['unresolved']))
 write.csv(do.call(rbind,out),'validation/data/tau_bias_audit_20261104/production_timing.csv',row.names=FALSE)
 stopifnot(z['unresolved']==0)
}
parallel::stopCluster(cl)
cat('8-worker/180-point linear estimate, hours:',sum(vapply(out,function(x)x$seconds,0.))*3*10000/32*2/8/3600,'\n')
