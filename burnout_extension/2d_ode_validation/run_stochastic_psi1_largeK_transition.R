source('psi_validation/R/stochastic_psi.R')
stopifnot(requireNamespace('data.table',quietly=TRUE),
  requireNamespace('adaptivetau',quietly=TRUE))

dir.create('2d_ode_validation/data',FALSE,TRUE)
psi<-1;theta<-.5;tau_epsilon<-.01;batch_size<-1000L
base_attempts<-1000L;precision_tiers<-c(1000L,3000L,10000L);ci_width_threshold<-.05
grid<-expand.grid(rho=c(.01,.02,.05,.10),K=c(1e6,1e7,1e8,1e9),
  q=seq(1.5,6,by=.25),KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE)
# For psi=1, theta=1/2 and high R0, K*rho^2/R0 is the leading trough-size
# scale.  q rescales it by rho*sqrt(K), so R0=rho*sqrt(K)/q follows the
# moving transition instead of holding R0 fixed as K grows.
grid$R0<-grid$rho*sqrt(grid$K)/grid$q
grid<-grid[grid$R0>1.05,]
grid<-grid[order(grid$rho,grid$K,grid$q),]
grid$theta<-theta;grid$psi<-psi;grid$point_id<-seq_len(nrow(grid))
checkpoint<-'2d_ode_validation/data/psi1_largeK_transition_tau_checkpoint.rds'
outfile<-'2d_ode_validation/data/psi1_largeK_transition_tau_results.csv'
version<-'psi1_largeK_transition_q15_6_by025_adaptivetau_eps001_v1'

new_state<-function()list(version=version,grid=grid,counts=transform(grid,
  method='adaptive tau',tau_epsilon=tau_epsilon,attempts=0L,
  established=0L,persistent=0L,unresolved=0L,target_attempts=base_attempts,
  batches=0L,elapsed_seconds=0))
state<-if(file.exists(checkpoint))readRDS(checkpoint)else new_state()
stopifnot(identical(state$version,version),identical(state$grid,grid))
redo<-which(state$counts$unresolved>0L)
if(length(redo)){
  fresh<-new_state()$counts
  state$counts[redo,names(fresh)]<-fresh[redo,names(fresh)]
}

save_state<-function(){
  tmp<-paste0(checkpoint,'.tmp');saveRDS(state,tmp,compress=FALSE)
  if(!file.copy(tmp,checkpoint,TRUE))stop('checkpoint replace failed')
  file.remove(tmp)
  z<-data.table::as.data.table(state$counts)
  z[,`:=`(P_unconditional=persistent/attempts,uncond_low=NA_real_,
    uncond_high=NA_real_,P_conditional=NA_real_,cond_low=NA_real_,cond_high=NA_real_)]
  for(j in which(z$attempts>0L)){
    u<-wilson(z$persistent[j],z$attempts[j]);z[j,`:=`(uncond_low=u[1],uncond_high=u[2])]
    if(z$established[j]>0L){cc<-wilson(z$persistent[j],z$established[j])
      z[j,`:=`(P_conditional=persistent/established,cond_low=cc[1],cond_high=cc[2])]}
  }
  data.table::fwrite(z,outfile)
}

cl<-parallel::makeCluster(16L);on.exit(parallel::stopCluster(cl),add=TRUE)
parallel::clusterExport(cl,'one_adaptive_tau_psi_chunked',envir=environment())
for(i in seq_len(nrow(grid)))repeat{
  x<-state$counts[i,]
  if(x$attempts>=x$target_attempts){
    wu<-diff(wilson(x$persistent,x$attempts))
    wc<-if(x$established>0)diff(wilson(x$persistent,x$established))else Inf
    if(max(wu,wc)>ci_width_threshold){
      higher<-precision_tiers[precision_tiers>x$attempts]
      if(length(higher)){state$counts$target_attempts[i]<-higher[1L];save_state();next}
    }
    break
  }
  n<-min(batch_size,x$target_attempts-x$attempts)
  seed<-1980000000L+i*100L+x$batches+1L;t0<-proc.time()[['elapsed']]
  got<-simulate_counts_batch_tau_psi_chunked_lb(cl,x$R0,x$rho,theta,psi,x$K,n,
    seed,tau_epsilon=tau_epsilon,chunk_time=50,tasks_per_worker=4L)
  state$counts$attempts[i]<-x$attempts+got['attempts']
  state$counts$established[i]<-x$established+got['established']
  state$counts$persistent[i]<-x$persistent+got['persistent']
  state$counts$unresolved[i]<-x$unresolved+got['unresolved']
  state$counts$batches[i]<-x$batches+1L
  state$counts$elapsed_seconds[i]<-x$elapsed_seconds+proc.time()[['elapsed']]-t0
  save_state();cat(sprintf('%d/%d q=%.2f R0=%.5g tau %d/%d\n',i,nrow(grid),x$q,x$R0,
    state$counts$attempts[i],state$counts$target_attempts[i]))
  if(got['unresolved']>0)stop('unresolved at point ',i)
}
save_state();cat('psi=1 large-K transition scan complete\n')
