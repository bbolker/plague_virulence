source('psi_validation/R/stochastic_psi.R')
stopifnot(requireNamespace('data.table',quietly=TRUE))
stopifnot(requireNamespace('Rcpp',quietly=TRUE))
Rcpp::sourceCpp('psi_validation/R/ctmc_psi.cpp')

dir.create('2d_ode_validation/data',FALSE,TRUE)
psi<-1;batch_size<-500L;base_attempts<-3000L;expanded_attempts<-10000L
grid<-expand.grid(rho=c(.01,.02,.05,.10),theta=c(0,.5,1),
  K=c(1000,3000,10000,30000),R0=1+c(.05,.075,.10,.15,.20,.30,.50,.75,1,1.5,2,3,5))
grid$scan<-'exact';grid<-grid[order(grid$rho,grid$theta,grid$K,grid$R0),]
grid$psi<-psi;grid$point_id<-seq_len(nrow(grid))
checkpoint<-'2d_ode_validation/data/psi1_scan_checkpoint.rds'
outfile<-'2d_ode_validation/data/psi1_stochastic_results.csv'
version<-'psi1_dynamic_g_v1'

new_state<-function()list(version=version,grid=grid,counts=transform(grid,
  attempts=0L,established=0L,persistent=0L,unresolved=0L,
  target_attempts=base_attempts,batches=0L,elapsed_seconds=0))
state<-if(file.exists(checkpoint))readRDS(checkpoint)else new_state()
stopifnot(identical(state$version,version),identical(state$grid,grid))

save_state<-function(){
  tmp<-paste0(checkpoint,'.tmp')
  saveRDS(state,tmp,compress=FALSE);file.copy(tmp,checkpoint,TRUE);file.remove(tmp)
  z<-data.table::as.data.table(state$counts)
  z[,`:=`(P_unconditional=persistent/attempts,uncond_low=NA_real_,
    uncond_high=NA_real_,P_conditional=NA_real_,cond_low=NA_real_,cond_high=NA_real_)]
  for(j in which(z$attempts>0)){
    u<-wilson(z$persistent[j],z$attempts[j])
    z[j,`:=`(uncond_low=u[1],uncond_high=u[2])]
    if(z$established[j]>0){
      cc<-wilson(z$persistent[j],z$established[j])
      z[j,`:=`(P_conditional=persistent/established,cond_low=cc[1],cond_high=cc[2])]
    }
  }
  data.table::fwrite(z,outfile)
}

for(i in seq_len(nrow(grid)))repeat{
  x<-state$counts[i,]
  if(x$attempts>=x$target_attempts){
    if(x$attempts==base_attempts){
      wu<-diff(wilson(x$persistent,x$attempts))
      wc<-if(x$established>0)diff(wilson(x$persistent,x$established))else Inf
      if(max(wu,wc)>.05){state$counts$target_attempts[i]<-expanded_attempts;save_state();next}
    }
    break
  }
  n<-min(batch_size,x$target_attempts-x$attempts)
  seed<-1960000000L+i*100L+x$batches+1L
  t0<-proc.time()[['elapsed']]
  got<-simulate_counts_cpp(x$R0,x$rho,x$theta,psi,x$K,n,seed)
  state$counts$attempts[i]<-x$attempts+got['attempts']
  state$counts$established[i]<-x$established+got['established']
  state$counts$persistent[i]<-x$persistent+got['persistent']
  state$counts$unresolved[i]<-x$unresolved+got['unresolved']
  state$counts$batches[i]<-x$batches+1L
  state$counts$elapsed_seconds[i]<-x$elapsed_seconds+proc.time()[['elapsed']]-t0
  save_state()
  cat(sprintf('%d/%d exact %d/%d\n',i,nrow(grid),state$counts$attempts[i],state$counts$target_attempts[i]))
  if(got['unresolved']>0)stop('unresolved at point ',i)
}
save_state();cat('psi=1 scan complete\n')
