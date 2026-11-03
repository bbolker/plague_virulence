source('validation/R/theory.R');source('validation/R/extrema.R');source('validation/R/fold_geometry.R');library(data.table)
# Research diagnostic: test the narrow theta/rho window exposed by the fold
# parameterization, rather than repeating the broad logarithmic parameter scan.
bs<-exp(seq(log(1e-12),log(100),length.out=100000));sm<-ifelse(bs<1e-3,bs/2+bs^2/6+bs^3/24+bs^4/120+bs^5/720,expm1(bs)/bs-1)
lookup<-approxfun(log(sm),log(bs),rule=2)
FofS<-function(s){B<-exp(lookup(log(s-1)));small<-B<1e-3;D<-exp(B)-s;D[small]<-B[small]/2+B[small]^2/3+B[small]^3/8+B[small]^4/30;s*D}
rr<-seq(1.1,25,length.out=9000);ths<-seq(.53,.60,.005);rhos<-seq(.24,.40,by=.0001)
cache<-'validation/data/narrow_fold_constituents.rds'
if(file.exists(cache))base<-readRDS(cache)else{
 base<-lapply(ths,function(th){z<-rbindlist(lapply(rr,function(r){cd<-C_log_derivatives(r,th);ad<-action_barrier_derivatives(r,th);q<-r*(r-1);qp<-2*r-1;c0<-cd['d1']+.5/(r-1)-th/(2*r);c1<-cd['d2']-.5/(r-1)^2+th/(2*r^2);data.table(R0=r,U=unname(-q*c0),V=unname(q*ad['d1']),U1=unname(-q*(qp*c0+q*c1)),V1=unname(q*(qp*ad['d1']+q*ad['d2'])),logC=unname(cd['logC']),A=unname(ad['value']))}));z[,theta:=th];cat('constituents',th,'\n');flush.console();z});saveRDS(base,cache)
}
out<-list();jj<-0;maxima<-list();kk<-0
for(a in base){th<-a$theta[1]
 for(rho in rhos){s<-a$U+a$V/rho;st<-a$U1+a$V1/rho;keep<-which(s>1);z<-rep(NA_real_,length(s));z[keep]<-FofS(s[keep])+st[keep]
  ix<-which(head(z,-1)*tail(z,-1)<0);n<-length(ix)
  kk<-kk+1;maxima[[kk]]<-data.table(theta=th,rho=rho,nturns=n)
  if(n){for(i in ix){r<-a$R0[i]-z[i]*(a$R0[i+1]-a$R0[i])/(z[i+1]-z[i]);sv<-approx(a$R0,s,xout=r)$y;B<-exp(lookup(log(sv-1)));lc<-approx(a$R0,a$logC,xout=r)$y;Av<-approx(a$R0,a$A,xout=r)$y;logb<-lc+.5*(log(rho)+log(r-1)-th*log(r)-log(2*pi))-Av/rho
   jj<-jj+1;out[[jj]]<-data.table(theta=th,rho=rho,R0=r,logK=log(B)-logb,direction=if(z[i]>0)'maximum'else 'minimum')
  }}
 };cat('classified',th,'\n');flush.console()
}
m<-rbindlist(maxima);d<-rbindlist(out);fwrite(m,'validation/data/narrow_fold_geometry_counts.csv');fwrite(d,'validation/data/narrow_fold_geometry_turns.csv');print(m[,.(maxturns=max(nturns)),by=theta]);print(m[nturns>2][1:10])
