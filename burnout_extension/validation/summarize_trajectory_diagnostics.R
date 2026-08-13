fs <- list.files('validation/data/trajectory_diagnostics',pattern='[.]rds$',full.names=TRUE)
for(f in fs){
  z<-readRDS(f);s<-z$summaries
  cat(basename(f),'\n');print(table(s$outcome))
  cat('first_trough=',sum(s$first_trough),' next_peak=',sum(s$outcome=='next_peak'),
      ' post_trough_extinction=',sum(s$outcome=='post_trough_extinction'),'\n')
  if(any(s$first_trough))print(summary(s$I_trough[s$first_trough]))
}
