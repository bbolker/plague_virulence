source('validation/R/stochastic.R'); source('validation/R/scan.R')
pts<-data.frame(
  R0=c(3.041169,3.041169,5.286084,3.04116887138897,3.9578069820069),
  rho=c(.04,.06,.10,.02,.02),theta=c(.25,.25,.50,0,0),
  K=c(3000,1000,3000,30000,30000))
out<-vector('list',nrow(pts))
for(i in seq_len(nrow(pts))) {
  cat('stochastic point',i,'/',nrow(pts),'\n')
  s<-do.call(simulate_point_unconditional,c(as.list(pts[i,]),list(n_attempts=10000,seed=20260812+i)))
  p<-do.call(run_point,as.list(pts[i,]))
  p_not_fizzle<-1-1/pts$R0[i]
  s$P1_uncond_K_ODE<-p_not_fizzle*p$P1_ref
  s$P1_uncond_K_first<-p_not_fizzle*p$P1_K_first
  s$P1_uncond_K_second<-p_not_fizzle*p$P1_K_second
  s$P1_uncond_L_second<-p_not_fizzle*p$P1_L_second
  out[[i]]<-s
}
data.table::fwrite(data.table::rbindlist(out),'validation/data/stochastic_results.csv')

# Complete event histories for a diagnostic sample from every parameter point.
dir.create('validation/data/trajectory_diagnostics',FALSE,TRUE)
old<-list.files('validation/data/trajectory_diagnostics',full.names=TRUE)
if(length(old)) unlink(old)
for(i in seq_len(nrow(pts))){
  set.seed(20261813+i);tr<-vector('list',100L);sm<-vector('list',100L)
  for(j in seq_len(100L)){
    z<-do.call(one_ctmc_trace,as.list(pts[i,]));tr[[j]]<-z$trajectory
    z$summary$trajectory_id<-j;sm[[j]]<-z$summary
  }
  tag<-sprintf('R0_%0.6f_rho_%0.2f_theta_%g_K_%d',pts$R0[i],pts$rho[i],pts$theta[i],pts$K[i])
  saveRDS(list(parameters=pts[i,],summaries=data.table::rbindlist(sm),trajectories=tr),
          file.path('validation/data/trajectory_diagnostics',paste0(tag,'.rds')),compress='xz')
}
