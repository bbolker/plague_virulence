source(file.path("fadeout","logistic_burnout","logistic_burnout_functions.R"))
out <- file.path("fadeout","logistic_burnout","outputs")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
cases <- data.frame(R0=c(1.5,2,2.5,3,5),r=c(.05,.125,.125,.2,.3))
rows <- list()
for(i in seq_len(nrow(cases))) {
  p <- cases[i,]
  z <- logistic_multitrough_probabilities(p$R0,p$r,10000,max_troughs=5,dt=.02)
  zf <- logistic_multitrough_probabilities(p$R0,p$r,10000,max_troughs=5,dt=.01)
  old <- logistic_burnout_probability(p$R0,p$r,10000,dt=.02,
                                      initial_tmax=100,maximum_tmax=5000)
  if(!nrow(z$trough_table)) next
  tt <- z$trough_table
  fine <- zf$trough_table
  n <- min(nrow(tt),nrow(fine))
  resolution <- rep(NA_real_,nrow(tt))
  resolution[seq_len(n)] <- abs(tt$P_persistence[seq_len(n)]-
                                  fine$P_persistence[seq_len(n)])
  identity_error <- abs(sum(tt$burnout_at_this_trough)+
                          tail(tt$cumulative_persistence,1)-1)
  rows[[i]] <- data.frame(
    R0=p$R0,r=p$r,trough_index=tt$trough_index,
    first_x_compatible=abs(tt$x_in[1]-old$x_in)<2e-7,
    first_q_compatible=abs(tt$q_lineage[1]-old$q1)<2e-7,
    first_Q_compatible=abs(tt$Q_burnout[1]-old$P_burnout)<2e-7,
    probabilities_valid=tt$q_lineage>=0&tt$q_lineage<=1&
      tt$Q_burnout>=0&tt$Q_burnout<=1&tt$P_persistence>=0&
      tt$P_persistence<=1&tt$cumulative_persistence>=0&
      tt$cumulative_persistence<=1&tt$burnout_at_this_trough>=0&
      tt$burnout_at_this_trough<=1,
    identity_error=identity_error,
    event_ordered=tt$t_peak<tt$t_in&
      c(TRUE,diff(tt$t_peak)>0)&c(TRUE,diff(tt$t_in)>0),
    dt_P_difference=resolution,
    x_nondecreasing=c(TRUE,diff(tt$x_in)>=-1e-8),
    P_nondecreasing=c(TRUE,diff(tt$P_persistence)>=-1e-8),
    Q_nonincreasing=c(TRUE,diff(tt$Q_burnout)<=1e-8),
    termination_status=z$termination_status
  )
}
validation <- do.call(rbind,rows)
write.csv(validation,file.path(out,"multitrough_validation_results.csv"),
          row.names=FALSE)
hard_ok <- with(validation,first_x_compatible&first_q_compatible&
  first_Q_compatible&probabilities_valid&identity_error<1e-8&
  event_ordered&(is.na(dt_P_difference)|dt_P_difference<2e-4))
if(!all(hard_ok)) stop("Hard multi-trough validation failure")
cat("Multi-trough validation rows: ",nrow(validation),"; all hard checks passed\n",
    sep="")
cat("Qualitative violations: x=",sum(!validation$x_nondecreasing),
    ", P=",sum(!validation$P_nondecreasing),
    ", Q=",sum(!validation$Q_nonincreasing),"\n",sep="")
