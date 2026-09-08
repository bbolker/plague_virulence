source('validation/R/stochastic.R')

one_ctmc_psi <- function(R0,rho,theta,psi,K,max_events=2e7) {
  S<-K-1L;I<-1L; first_peak<-FALSE;first_trough<-FALSE
  s_drift_positive<-FALSE;s_turned<-FALSE
  growth <- function(S,I) R0*(S/K)*((S+I)/K)^(-psi)-1
  for(ev in seq_len(max_events)) {
    Sold<-S;Iold<-I; gold<-growth(S,I)
    inf<-R0*S*I/K*(K/(S+I))^psi;rem<-I
    rec<-if(S<K) K*rho*(1-S/K)*(S/K)^theta else 0
    tot<-inf+rem+rec;if(tot<=0)break
    u<-runif(1)*tot
    if(u<inf){S<-S-1L;I<-I+1L}else if(u<inf+rem)I<-I-1L else S<-S+1L
    gnew<-if(I>0) growth(S,I) else NA_real_
    if(!first_peak&&gold>0&&is.finite(gnew)&&gnew<=0)first_peak<-TRUE
    if(first_peak&&!first_trough&&gold<0&&is.finite(gnew)&&gnew>=0)first_trough<-TRUE
    if(first_trough){
      drift<-K*rho*(1-S/K)*(S/K)^theta-
        R0*S*I/K*(K/(S+I))^psi
      if(drift>0)s_drift_positive<-TRUE
      if(s_drift_positive&&drift<0)s_turned<-TRUE
    }
    if(I==0)return(c(established=first_peak,persist=FALSE))
    if(first_trough&&s_turned&&gold>0&&gnew<=0)
      return(c(established=TRUE,persist=TRUE))
  }
  c(established=first_peak,persist=NA)
}

simulate_counts_batch_psi <- function(cl,R0,rho,theta,psi,K,n_attempts,seed) {
  nw<-length(cl);chunks<-rep(n_attempts%/%nw,nw)
  if(n_attempts%%nw)chunks[seq_len(n_attempts%%nw)]<-chunks[seq_len(n_attempts%%nw)]+1L
  parallel::clusterSetRNGStream(cl,iseed=seed)
  ans<-parallel::parLapply(cl,chunks,function(n,R0,rho,theta,psi,K){
    out<-c(attempts=n,established=0L,persistent=0L,unresolved=0L)
    for(j in seq_len(n)){
      z<-one_ctmc_psi(R0,rho,theta,psi,K)
      if(isTRUE(z['established']))out['established']<-out['established']+1L
      if(is.na(z['persist']))out['unresolved']<-out['unresolved']+1L
      else if(isTRUE(z['persist']))out['persistent']<-out['persistent']+1L
    };out
  },R0=R0,rho=rho,theta=theta,psi=psi,K=K)
  colSums(do.call(rbind,ans))
}

one_adaptive_tau_psi <- function(R0,rho,theta,psi,K,tf=max(500,20/rho),
                                 tau_epsilon=.01) {
  if(!requireNamespace('adaptivetau',quietly=TRUE))stop('adaptivetau required')
  transitions<-list(c(S=-1,I=1),c(I=-1),c(S=1))
  rates<-function(z,p,t){S<-z['S'];I<-z['I'];N<-S+I
    c(infection=if(I>0) p$R0*S*I/p$K*(p$K/N)^p$psi else 0,
      removal=I,recruitment=if(S<p$K)p$K*p$rho*(1-S/p$K)*(S/p$K)^p$theta else 0)}
  z<-adaptivetau::ssa.adaptivetau(c(S=K-1,I=1),transitions,rates,
    list(R0=R0,rho=rho,theta=theta,psi=psi,K=K),tf=tf,
    tl.params=list(epsilon=tau_epsilon,extraChecks=TRUE))
  S<-z[,'S'];I<-z[,'I'];N<-S+I
  gg<-R0*(S/K)*(K/N)^psi-1
  down<-which(gg[-length(gg)]>0 & gg[-1L]<=0 & I[-1L]>0)+1L
  if(!length(down))return(c(established=FALSE,persist=FALSE))
  peak<-down[1L]; after<-seq.int(peak+1L,length(S))
  up<-after[gg[after-1L]<0 & gg[after]>=0 & I[after]>0]
  extinct<-which(I==0)
  if(!length(up)){
    if(length(extinct))return(c(established=TRUE,persist=FALSE))
    return(c(established=TRUE,persist=NA))
  }
  trough<-up[1L]
  if(length(extinct)&&extinct[1L]>trough)return(c(established=TRUE,persist=FALSE))
  if(any(down>trough))return(c(established=TRUE,persist=TRUE))
  c(established=TRUE,persist=NA)
}

simulate_counts_batch_tau_psi <- function(cl,R0,rho,theta,psi,K,n_attempts,seed) {
  nw<-length(cl);chunks<-rep(n_attempts%/%nw,nw)
  if(n_attempts%%nw)chunks[seq_len(n_attempts%%nw)]<-chunks[seq_len(n_attempts%%nw)]+1L
  parallel::clusterSetRNGStream(cl,iseed=seed)
  ans<-parallel::parLapply(cl,chunks,function(n,R0,rho,theta,psi,K){
    out<-c(attempts=n,established=0L,persistent=0L,unresolved=0L)
    for(j in seq_len(n)){
      z<-one_adaptive_tau_psi(R0,rho,theta,psi,K)
      if(isTRUE(z['established']))out['established']<-out['established']+1L
      if(is.na(z['persist']))out['unresolved']<-out['unresolved']+1L
      else if(isTRUE(z['persist']))out['persistent']<-out['persistent']+1L
    };out
  },R0=R0,rho=rho,theta=theta,psi=psi,K=K)
  colSums(do.call(rbind,ans))
}
