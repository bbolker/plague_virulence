source('validation/R/theory.R');source('validation/R/extrema.R');library(data.table)
ans<-rbindlist(lapply(c(0,.25,.5,.75,.9,.95,.975),function(th){
 f<-function(r)unname(action_barrier_derivatives(r,th)['d1'])
 root<-uniroot(f,c(1.01,600),tol=1e-10)$root
 data.table(theta=th,action_maximum_R0=root)
}));print(ans);fwrite(ans,'validation/data/action_barrier_sign_boundaries.csv')
# Check the derived identity for rho^2*dE/drho using differences in rho.
rows<-rbindlist(lapply(c(1.2,2,3,6,10),function(r)rbindlist(lapply(c(0,.5,1),function(th){
 cd<-C_log_derivatives(r,th);ad<-action_barrier_derivatives(r,th);a<-ad['value'];ap<-ad['d1'];q<-r*(r-1);c0<-cd['d1']+1/(2*(r-1))-th/(2*r);K<-1e6
 E<-function(rv){lb<-log(K)+cd['logC']+.5*(log(rv)+log(r-1)-th*log(r)-log(2*pi))-a/rv;B<-exp(lb);if(B>50)return(NA_real_);H<-if(B<1e-4)1+B/2+B^2/6 else expm1(B)/B;c0-ap/rv+H/q}
 rbindlist(lapply(c(.001,.003,.01,.03),function(rv){lb<-log(K)+cd['logC']+.5*(log(rv)+log(r-1)-th*log(r)-log(2*pi))-a/rv;B<-exp(lb);if(B>40)return(NULL);D<-if(B<1e-4)B/2+B^2/3+B^3/8 else ((B-1)*exp(B)+1)/B
 exact<-ap+D*(rv/2+a)/q;h<-rv*1e-5;fd<-rv^2*(E(rv+h)-E(rv-h))/(2*h);data.table(R0=r,theta=th,rho=rv,identity=as.numeric(exact),finite_difference=as.numeric(fd),scaled_error=abs(exact-fd)/max(1,abs(exact)))
 }))
}))))
fwrite(rows,'validation/data/rho_structure_identity_checks.csv');print(max(rows$scaled_error,na.rm=TRUE));stopifnot(max(rows$scaled_error,na.rm=TRUE)<1e-5)
