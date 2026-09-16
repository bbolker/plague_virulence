source('2d_ode_validation/R/trough_kendall_psi.R')
stopifnot(requireNamespace('data.table',quietly=TRUE),requireNamespace('adaptivetau',quietly=TRUE))
library(data.table)

R0<-4;rho<-.01;theta<-.5;psi<-1;K<-1e6
z<-trough_kendall_2d(R0,rho,theta,psi,K);stopifnot(z$status=='OK')
T<-z$t_plus-z$t_minus
S0<-as.integer(round(K*z$x_minus));I0<-as.integer(round(K*z$y_minus))

one_peak_full<-function(seed,S0,I0,T,R0,rho,theta,psi,K,epsilon=.01,chunk=5){
  set.seed(seed)
  transitions<-list(c(S=-1,I=1),c(I=-1),c(S=1))
  rates<-function(v,p,t){S<-v['S'];I<-v['I'];N<-S+I
    c(infection=if(I>0)p$R0*S*I/p$K*(p$K/N)^p$psi else 0,
      removal=I,recruitment=if(S<p$K)p$K*p$rho*(1-S/p$K)*(S/p$K)^p$theta else 0)}
  state<-c(S=S0,I=I0);elapsed<-0
  while(elapsed<T&&state['I']>0){
    h<-min(chunk,T-elapsed)
    a<-adaptivetau::ssa.adaptivetau(state,transitions,rates,
      list(R0=R0,rho=rho,theta=theta,psi=psi,K=K),tf=h,
      tl.params=list(epsilon=epsilon,extraChecks=TRUE))
    state<-c(S=tail(a[,'S'],1),I=tail(a[,'I'],1));elapsed<-elapsed+h
  }
  state['I']>0
}

nrep<-5000L
cl<-parallel::makeCluster(16L);on.exit(parallel::stopCluster(cl),add=TRUE)
parallel::clusterEvalQ(cl,library(adaptivetau))
parallel::clusterExport(cl,'one_peak_full',envir=environment())
ids<-seq_len(nrep)
alive<-unlist(parallel::parLapplyLB(cl,ids,function(i,S0,I0,T,R0,rho,theta,psi,K)
  one_peak_full(73000000L+i,S0,I0,T,R0,rho,theta,psi,K),S0=S0,I0=I0,T=T,
  R0=R0,rho=rho,theta=theta,psi=psi,K=K))
wilson<-function(x,n,zq=qnorm(.975)){den<-1+zq^2/n;ctr<-(x/n+zq^2/(2*n))/den
  hw<-zq*sqrt((x/n)*(1-x/n)/n+zq^2/(4*n^2))/den;c(low=ctr-hw,high=ctr+hw)}
w<-wilson(sum(alive),nrep)
Pk_round<- -expm1(I0*log1p(-z$single_lineage_u))
out<-data.table(method=c('Full CTMC from deterministic first peak','Kendall from deterministic first peak'),
  estimate=c(mean(alive),Pk_round),low=c(w[1],NA),high=c(w[2],NA),
  n=c(nrep,NA_integer_),R0=R0,rho=rho,theta=theta,psi=psi,K=K,
  S0=S0,I0=I0,T=T,x0=z$x_minus,y0=z$y_minus,y_over_x=z$y_minus/z$x_minus,
  feedback=psi*z$y_minus/(z$x_minus+z$y_minus))
fwrite(out,'2d_ode_validation/data/psi1_peak_handoff_isolation.csv')
print(out)
