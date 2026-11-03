# Corrected fixed-design high-R scan. Original raw scan stays untouched.
source('validation/R/stochastic.R');source('validation/R/native_simulation.R')
stopifnot(requireNamespace('data.table',quietly=TRUE))
checkpoint<-'validation/data/fig1_high_R_guarded_checkpoint.rds'
output<-'validation/data/fig1_high_R_guarded_results.csv'
grid<-expand.grid(rho=.01,theta=c(0,.25,.5,.75,1),K=c(1e6,1e7,1e8,1e9),
                  R0=c(8,12,16,24,32,48,64,96,128),KEEP.OUT.ATTRS=FALSE)
grid<-grid[order(grid$theta,grid$K,grid$R0),];rownames(grid)<-NULL
grid$point_id<-seq_len(nrow(grid))
fingerprint<-tools::md5sum(c('validation/src/guarded_first_trough.c','validation/R/native_simulation.R'))
version<-'boundary_guard_v1_eps001_budget001_small100_fixed10000'
new_state<-function()list(version=version,fingerprint=fingerprint,grid=grid,
 seed_base=2026140000L,batch_size=500L,workers=8L,created=as.character(Sys.time()),
 counts=transform(grid,method='boundary-guarded adaptive tau + SSA',tau_epsilon=.01,
  event_budget=.01,small_count=100,tf=500,attempts=0L,established=0L,persistent=0L,
  unresolved=0L,ever_zero=0L,trough=0L,events=0,tau_steps=0,exact_steps=0,rejected=0,
  target_attempts=10000L,completed_batches=0L,elapsed_seconds=0))
state<-if(file.exists(checkpoint))readRDS(checkpoint) else new_state()
stopifnot(identical(state$version,version),identical(state$fingerprint,fingerprint),identical(state$grid,grid))
save_state<-function(){
 state$updated<<-as.character(Sys.time())
 tmp<-paste0(checkpoint,'.tmp');saveRDS(state,tmp,compress=FALSE)
 stopifnot(file.rename(tmp,checkpoint))
 z<-data.table::as.data.table(state$counts)
 z[,`:=`(P_unconditional=persistent/attempts,P_conditional=persistent/established,
          uncond_low=NA_real_,uncond_high=NA_real_,cond_low=NA_real_,cond_high=NA_real_,
          sample_source='boundary_guarded_high_R_fixed10000')]
 for(i in which(z$attempts>0)){
  a<-wilson(z$persistent[i],z$attempts[i]);b<-wilson(z$persistent[i],z$established[i])
  z[i,`:=`(uncond_low=max(0,a[1]),uncond_high=min(1,a[2]),cond_low=max(0,b[1]),cond_high=min(1,b[2]))]
 }
 tmp<-paste0(output,'.tmp');data.table::fwrite(z,tmp);stopifnot(file.rename(tmp,output))
}
dll<-compile_native_simulation('validation/src/guarded_first_trough.c')
cl<-parallel::makeCluster(8L);parallel::clusterCall(cl,function(p)dyn.load(p),dll)
save_state()
fields<-c('attempts','established','persistent','unresolved','ever_zero','trough','events','tau_steps','exact_steps','rejected')
for(i in seq_len(nrow(grid)))repeat{
 row<-state$counts[i,];if(row$attempts>=row$target_attempts)break
 row$method<-'native_guard';row$epsilon<-row$tau_epsilon
 batch<-row$completed_batches+1L;t<-proc.time()[3]
 got<-native_batch(cl,row,500L,state$seed_base+i*100L+batch)
 for(nm in fields)state$counts[[nm]][i]<-state$counts[[nm]][i]+got[nm]
 state$counts$completed_batches[i]<-batch
 state$counts$elapsed_seconds[i]<-row$elapsed_seconds+proc.time()[3]-t
 save_state();cat('point',i,'/180:',state$counts$attempts[i],'/10000\n');flush.console()
 if(got['unresolved']>0)stop('Unresolved trajectories; no figure replacement permitted.')
}
stopifnot(all(state$counts$attempts==10000L),all(state$counts$unresolved==0L))
state$completed<-as.character(Sys.time());save_state();parallel::stopCluster(cl)
cat('Corrected high-R scan complete.\n')
