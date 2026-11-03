source('validation/R/theory.R');source('validation/R/extrema.R');library(data.table)
logH<-function(B){if(B<1e-4)log1p(B/2+B^2/6+B^3/24)else if(B<700)log(expm1(B))-log(B)else B-log(B)}
Bc<-function(s){stopifnot(s>1); uniroot(function(v)logH(v)-log(s),c(1e-12,2*log(s)+10),tol=1e-12)$root}
stableP<-function(r,rv,th,k){cd<-C_log_derivatives(r,th);ad<-action_barrier_derivatives(r,th);lb<-log(k)+cd['logC']+.5*(log(rv)+log(r-1)-th*log(r)-log(2*pi))-ad['value']/rv;(1-1/r)*(-expm1(-exp(min(lb,700))))}
grid<-CJ(R0=c(1.1,1.2,1.5,2,3,5,10,50,300),rho=c(.002,.01,.03),theta=c(0,.25,.5,.75,1))
ans<-rbindlist(lapply(seq_len(nrow(grid)),function(i){
 r<-grid$R0[i];rv<-grid$rho[i];th<-grid$theta[i];q<-extremum_quantities(r,rv,th,1);s<--r*(r-1)*q['L'];if(s<=1)return(NULL)
 bc<-Bc(s);lk<-log(bc)-q['logB'];
 # Check sign reversal against direct symmetric differences of the original BI probability.
 rbindlist(lapply(c(.5,2),function(mult){
  if(lk+log(mult)>650)return(NULL);k<-exp(lk)*mult;h<-min(r*1e-5,(r-1)*1e-3)
  fd<-(stableP(r+h,rv,th,k)-stableP(r-h,rv,th,k))/(2*h)
  data.table(R0=r,rho=rv,theta=th,logKcritical=as.numeric(lk),factor=mult,finite_difference=fd,expected_sign=if(mult<1)-1 else 1)
 }))
}))
fwrite(ans,'validation/data/analytic_monotonicity_threshold_checks.csv');use<-ans[abs(finite_difference)>1e-9];print(use[,.(checks=.N,sign_mismatches=sum(sign(finite_difference)!=expected_sign))]);stopifnot(all(sign(use$finite_difference)==use$expected_sign))
