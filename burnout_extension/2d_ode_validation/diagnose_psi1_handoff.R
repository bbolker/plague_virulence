source('2d_ode_validation/R/trough_kendall_psi.R')
stopifnot(requireNamespace('data.table',quietly=TRUE),requireNamespace('deSolve',quietly=TRUE))
library(data.table)

R0<-4;rho<-.01;theta<-.5;psi<-1;K<-1e6
z<-trough_kendall_2d(R0,rho,theta,psi,K)
stopifnot(z$status=='OK')
T<-z$t_plus-z$t_t
S0<-as.integer(round(K*z$x_t));I0<-as.integer(round(K*z$y_t))

rhs<-function(t,s,p){x<-s[1];y<-s[2];n<-x+y
  list(c(rho*(1-x)*x^theta-R0*x*y/n^psi,(R0*x/n^psi-1)*y))}
tt<-seq(0,T,length.out=10001)
path<-as.data.table(deSolve::ode(c(x=z$x_t,y=z$y_t),tt,rhs,NULL,
  method='lsoda',rtol=1e-11,atol=c(1e-13,1e-13),maxsteps=1e6))
path[,b:=R0*x/(x+y)^psi]
bfun<-approxfun(path$time,path$b,rule=2)
bmax<-max(path$b)*1.000001
A<-sum(diff(path$time)*(head(exp(-(log(path$y)-log(path$y[1]))),-1)+
  tail(exp(-(log(path$y)-log(path$y[1]))),-1))/2)
u<-1/(1+A)
P_kendall_handoff<- -expm1(I0*log1p(-u))

one_full<-function(seed,S0,I0,T,R0,rho,theta,psi,K){
  set.seed(seed);S<-S0;I<-I0;t<-0
  while(t<T&&I>0){
    N<-S+I
    inf<-if(I>0)R0*S*I/K*(K/N)^psi else 0
    rem<-I
    rec<-if(S<K)K*rho*(1-S/K)*(S/K)^theta else 0
    tot<-inf+rem+rec
    if(!is.finite(tot)||tot<=0)break
    dt<-rexp(1,tot);if(t+dt>T)break;t<-t+dt
    v<-runif(1)*tot
    if(v<inf){S<-S-1L;I<-I+1L}else if(v<inf+rem){I<-I-1L}else S<-S+1L
  }
  I>0
}

# Exact thinning for the time-inhomogeneous linear birth-death process with
# per-capita birth b(t), death 1, and deterministic prescribed background.
one_frozen<-function(seed,I0,T,b_time,b_value,bmax){
  set.seed(seed);I<-I0;t<-0;bfun<-approxfun(b_time,b_value,rule=2)
  while(t<T&&I>0){
    tc<-t+rexp(1,(bmax+1)*I);if(tc>T)break
    b<-bfun(tc)
    if(runif(1)<(b+1)/(bmax+1)){
      if(runif(1)<b/(b+1))I<-I+1L else I<-I-1L
    }
    t<-tc
  }
  I>0
}

nrep<-20000L
cl<-parallel::makeCluster(16L);on.exit(parallel::stopCluster(cl),add=TRUE)
parallel::clusterExport(cl,c('one_full','one_frozen'),envir=environment())
ids<-seq_len(nrep)
full<-unlist(parallel::parLapplyLB(cl,ids,function(i,S0,I0,T,R0,rho,theta,psi,K)
  one_full(71000000L+i,S0,I0,T,R0,rho,theta,psi,K),S0=S0,I0=I0,T=T,
  R0=R0,rho=rho,theta=theta,psi=psi,K=K))
frozen<-unlist(parallel::parLapplyLB(cl,ids,function(i,I0,T,b_time,b_value,bmax)
  one_frozen(72000000L+i,I0,T,b_time,b_value,bmax),I0=I0,T=T,
  b_time=path$time,b_value=path$b,bmax=bmax))

wilson<-function(x,n,z=qnorm(.975)){den<-1+z^2/n;ctr<-(x/n+z^2/(2*n))/den
  hw<-z*sqrt((x/n)*(1-x/n)/n+z^2/(4*n^2))/den;c(low=ctr-hw,high=ctr+hw)}
out<-rbindlist(list(
  data.table(method='Full CTMC from deterministic trough',success=sum(full),n=nrep,
    estimate=mean(full)),
  data.table(method='Frozen-background branching simulation',success=sum(frozen),n=nrep,
    estimate=mean(frozen)),
  data.table(method='Kendall formula on frozen background',success=NA_integer_,n=NA_integer_,
    estimate=P_kendall_handoff)))
out[,`:=`(low=NA_real_,high=NA_real_)]
for(i in which(!is.na(out$n))){w<-wilson(out$success[i],out$n[i]);out[i,`:=`(low=w[1],high=w[2])]}
out[,`:=`(R0=R0,rho=rho,theta=theta,psi=psi,K=K,S0=S0,I0=I0,T=T,
  x0=z$x_t,y0=z$y_t,y_over_x=z$y_t/z$x_t,
  feedback=psi*z$y_t/(z$x_t+z$y_t))]
fwrite(out,'2d_ode_validation/data/psi1_handoff_isolation.csv')
fwrite(path,'2d_ode_validation/data/psi1_handoff_frozen_path.csv')
print(out)
