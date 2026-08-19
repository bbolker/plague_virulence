wilson <- function(k,n,z=1.96) {
  p<-k/n; d<-1+z^2/n; c((p+z^2/(2*n)-z*sqrt(p*(1-p)/n+z^2/(4*n^2)))/d,
                         (p+z^2/(2*n)+z*sqrt(p*(1-p)/n+z^2/(4*n^2)))/d)
}

one_ctmc <- function(R0,rho,theta,K,max_events=2e7) {
  S<-K-1L; I<-1L; xs<-1/R0
  first_peak<-FALSE;first_trough<-FALSE;s_drift_positive<-FALSE;s_turned<-FALSE
  for(ev in seq_len(max_events)) {
    Sold <- S
    inf<-R0*S*I/K; rem<-I
    rec<-if(S<K) K*rho*(1-S/K)*(S/K)^theta else 0
    tot<-inf+rem+rec; if(tot<=0) break
    u<-runif(1)*tot
    if(u<inf) {S<-S-1L;I<-I+1L} else if(u<inf+rem) I<-I-1L else S<-S+1L
    if(!first_peak && Sold/K>xs && S/K<=xs && I>0) first_peak<-TRUE
    if(first_peak && !first_trough && Sold/K<xs && S/K>=xs && I>0) first_trough<-TRUE
    if(first_trough) {
      g<-K*rho*(1-S/K)*(S/K)^theta-R0*S*I/K
      if(g>0) s_drift_positive<-TRUE
      if(s_drift_positive && g<0) s_turned<-TRUE
    }
    if(I==0) return(c(established=first_peak,persist=FALSE))
    if(first_trough && s_turned && Sold/K>xs && S/K<=xs && I>0)
      return(c(established=TRUE,persist=TRUE))
  }
  c(established=first_peak,persist=NA)
}

# Diagnostic version: retain the complete event history and continue beyond the
# first trough until extinction or the next infectious peak (down-crossing x*).
one_ctmc_trace <- function(R0,rho,theta,K,max_events=2e7) {
  S<-K-1L; I<-1L; xs<-1/R0; t<-0; first_peak<-FALSE; first_trough<-FALSE
  s_drift_positive<-FALSE;s_turned<-FALSE
  tt<-numeric(10000L); SS<-integer(10000L); II<-integer(10000L); n<-1L
  tt[1]<-0;SS[1]<-S;II[1]<-I
  grow <- function() {length(tt)<<-length(tt)*2L;length(SS)<<-length(SS)*2L;length(II)<<-length(II)*2L}
  outcome<-'unresolved';I_trough<-NA_integer_;t_trough<-NA_real_
  for(ev in seq_len(max_events)) {
    Sold<-S; inf<-R0*S*I/K; rem<-I
    rec<-if(S<K) K*rho*(1-S/K)*(S/K)^theta else 0
    tot<-inf+rem+rec;if(tot<=0) break
    t<-t+rexp(1,tot);u<-runif(1)*tot
    if(u<inf){S<-S-1L;I<-I+1L}else if(u<inf+rem)I<-I-1L else S<-S+1L
    n<-n+1L;if(n>length(tt))grow();tt[n]<-t;SS[n]<-S;II[n]<-I
    if(!first_peak&&Sold/K>xs&&S/K<=xs&&I>0)first_peak<-TRUE
    if(first_peak&&!first_trough&&Sold/K<xs&&S/K>=xs&&I>0){
      first_trough<-TRUE;I_trough<-I;t_trough<-t
    }
    if(first_trough){
      g<-K*rho*(1-S/K)*(S/K)^theta-R0*S*I/K
      if(g>0)s_drift_positive<-TRUE
      if(s_drift_positive&&g<0)s_turned<-TRUE
    }
    if(I==0){outcome<-if(first_trough)'post_trough_extinction' else 'extinction';break}
    if(first_trough&&s_turned&&Sold/K>xs&&S/K<=xs){outcome<-'next_peak';break}
  }
  list(summary=data.frame(outcome=outcome,first_peak=first_peak,first_trough=first_trough,
       I_trough=I_trough,t_trough=t_trough,s_turned=s_turned,
       t_end=t,S_end=S,I_end=I,events=n-1L),
       trajectory=data.frame(t=tt[seq_len(n)],S=SS[seq_len(n)],I=II[seq_len(n)]))
}

simulate_point <- function(R0,rho,theta,K,n_established=1000,seed=1) {
  set.seed(seed); est<-pers<-attempt<-0L
  while(est<n_established && attempt<50*n_established) {
    z<-one_ctmc(R0,rho,theta,K); attempt<-attempt+1L
    if(isTRUE(z['established'])) {est<-est+1L;if(isTRUE(z['persist'])) pers<-pers+1L}
  }
  ci<-wilson(pers,est)
  data.frame(R0=R0,rho=rho,theta=theta,K=K,attempts=attempt,established=est,
             persistent=pers,P1_hat=pers/est,CI_low=ci[1],CI_high=ci[2])
}

# Unconditional persistence beyond the first trough, matching
# burnout::P1_prob(): early fizzles and post-outbreak burnout are failures.
simulate_point_unconditional <- function(R0,rho,theta,K,n_attempts=10000,seed=1,
                                         n_cores=min(8L,parallel::detectCores())) {
  n_cores <- max(1L,min(as.integer(n_cores),as.integer(n_attempts)))
  chunks <- rep(n_attempts %/% n_cores,n_cores)
  chunks[seq_len(n_attempts %% n_cores)] <- chunks[seq_len(n_attempts %% n_cores)]+1L
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl),add=TRUE)
  parallel::clusterExport(cl,'one_ctmc',envir=environment())
  parallel::clusterSetRNGStream(cl,iseed=seed)
  ans <- parallel::parLapply(cl,chunks,function(n,R0,rho,theta,K) {
    persistent <- established <- unresolved <- 0L
    for(attempt in seq_len(n)) {
      z <- one_ctmc(R0,rho,theta,K)
      if(isTRUE(z['established'])) established <- established+1L
      if(is.na(z['persist'])) unresolved <- unresolved+1L
      else if(isTRUE(z['persist'])) persistent <- persistent+1L
    }
    c(established=established,persistent=persistent,unresolved=unresolved)
  },R0=R0,rho=rho,theta=theta,K=K)
  totals <- colSums(do.call(rbind,ans))
  established <- totals['established']; persistent <- totals['persistent']
  unresolved <- totals['unresolved']
  if(unresolved>0L) stop('Unresolved Gillespie trajectories: ',unresolved)
  ci <- wilson(persistent,n_attempts)
  data.frame(R0=R0,rho=rho,theta=theta,K=K,attempts=n_attempts,
             established=established,persistent=persistent,
             P1_hat=persistent/n_attempts,CI_low=ci[1],CI_high=ci[2])
}

# Reusable-cluster batch for checkpointed scans.  Counts from independent
# batches add exactly, so a stopped scan can resume without rerunning completed
# batches.
simulate_counts_batch <- function(cl,R0,rho,theta,K,n_attempts,seed) {
  n_workers <- length(cl)
  chunks <- rep(n_attempts %/% n_workers,n_workers)
  chunks[seq_len(n_attempts %% n_workers)] <-
    chunks[seq_len(n_attempts %% n_workers)]+1L
  parallel::clusterSetRNGStream(cl,iseed=seed)
  ans <- parallel::parLapply(cl,chunks,function(n,R0,rho,theta,K) {
    persistent <- established <- unresolved <- 0L
    for(attempt in seq_len(n)) {
      z <- one_ctmc(R0,rho,theta,K)
      if(isTRUE(z['established'])) established <- established+1L
      if(is.na(z['persist'])) unresolved <- unresolved+1L
      else if(isTRUE(z['persist'])) persistent <- persistent+1L
    }
    c(attempts=n,established=established,persistent=persistent,
      unresolved=unresolved)
  },R0=R0,rho=rho,theta=theta,K=K)
  colSums(do.call(rbind,ans))
}
