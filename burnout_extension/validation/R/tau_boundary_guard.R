# Keep the expected number of events affecting a small population per leap
# small. This guards absorbing-boundary hitting probabilities, which are not
# controlled by the ordinary epsilon tolerance alone (its count bound floors
# at one). All transition rates and the outcome classifier remain unchanged.
tau_boundary_max_step <- function(z,p,t) {
  S<-unname(z['S']);I<-unname(z['I'])
  infection<-p$R0*S*I/p$K
  recruitment<-if(S<p$K)p$K*p$rho*(1-S/p$K)*(S/p$K)^p$theta else 0
  cap<-Inf
  if(p$theta>0 && S>0 && S<=p$small_count && infection+recruitment>0)
    cap<-min(cap,p$event_budget/(infection+recruitment))
  if(I>0 && I<=p$small_count && infection+I>0)
    cap<-min(cap,p$event_budget/(infection+I))
  cap
}

one_boundary_guarded_tau <- function(R0,rho,theta,K,tau_epsilon=.01,
                                    event_budget=.01,small_count=100,tf=100) {
  trans<-list(c(S=-1,I=1),c(I=-1),c(S=1))
  rates<-function(z,p,t){S<-z['S'];I<-z['I'];c(p$R0*S*I/p$K,I,
    if(S<p$K)p$K*p$rho*(1-S/p$K)*(S/p$K)^p$theta else 0)}
  params<-list(R0=R0,rho=rho,theta=theta,K=K,
               event_budget=event_budget,small_count=small_count)
  z<-adaptivetau::ssa.adaptivetau(c(S=K-1,I=1),trans,rates,params,tf=tf,
       maxTauFunc=tau_boundary_max_step,
       tl.params=list(epsilon=tau_epsilon,maxtau=0,extraChecks=TRUE))
  a<-classify_first_trough(z[,'S'],z[,'I'],R0,rho,theta,K)
  c(established=as.integer(a['established']),persistent=as.integer(isTRUE(a['persist'])),
    unresolved=as.integer(is.na(a['persist'])),ever_zero=as.integer(any(z[,'S']==0)))
}
