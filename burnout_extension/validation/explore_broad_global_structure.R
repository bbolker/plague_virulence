source('validation/R/theory.R');source('validation/R/extrema.R');library(data.table)
# Stable large-R evaluation retaining exp(-(1-theta)R), which is essential near theta=1.
large_C<-function(r,th){d<-1-th
 low<-function(s){ans<-numeric(length(s));small<-s<1e-4;ss<-s[small];ans[small]<-.5-d+(d*d/2-1/12)*ss-d^3*ss^2/6;ss<-s[!small];ans[!small]<-expm1(-d*ss)/ss+1/ss-1/expm1(ss);ans}
 ib<-integrate(low,0,1,rel.tol=1e-10)$value+integrate(function(v)exp(-d*exp(v))-ifelse(v<6,exp(v)/expm1(exp(v)),0),0,log(r),rel.tol=1e-10,subdivisions=400)$value
 T<-min(55,r-1);co<-integrate(function(t){s<-r-t;u<-exp(-t);den<--r*expm1(-t)-t;exp(-d*s)*((1-r*u)*(-expm1(-t))/den-1/s)},0,T,rel.tol=1e-9,abs.tol=1e-12)$value
 -log(r)+ib+co
}
large_A<-function(r,th){d<-1-th;n<-1:18
 if(d==0)return(r-log(r)-1+sum(r^-n/(n*(n+1))))
 t<-d*r
 lead<-if(t<40)exp(-t)*expm1(d*(r-log(r))-log1p(d))/d else exp(-d*log(r))/(d*(1+d))
 lead+exp(-d*log(r))*sum(r^-n/((n+d)*(n+d+1)))
}
large_dA<-function(r,th){d<-1-th;n<-0:18;exp(-d*r)-r^(-d-1)*sum(r^-n/(n+d+1))}
# Independent overlap checks against full quadrature/automatic derivatives.
checks<-rbindlist(lapply(c(50,100,300,600),function(r)rbindlist(lapply(c(0,.5,.975,.99,.999,1),function(th){cd<-C_log_derivatives(r,th);ad<-action_barrier_derivatives(r,th);h<-r*1e-4;dc<-(large_C(r+h,th)-large_C(r-h,th))/(2*h);data.table(R0=r,theta=th,logC_error=large_C(r,th)-cd['logC'],A_error=large_A(r,th)-ad['value'],r_dC_error=r*(dc-cd['d1']),dA_error=large_dA(r,th)-ad['d1'])}))))
fwrite(checks,'validation/data/uniform_largeR_crosschecks.csv');print(checks[,lapply(.SD,function(x)max(abs(x))),.SDcols=3:6])
ths<-c(seq(0,.95,.05),.975,.99,.995,.999,.9999,1)
rhos<-10^seq(-6,3,length.out=125);ks<-10^seq(0,15,by=1.5)
rr<-unique(c(1+exp(seq(log(1e-9),log(9),length.out=950)),exp(seq(log(10),log(600),length.out=200))))
cache<-'validation/data/broad_structure_base.rds'
if(file.exists(cache))bas<-readRDS(cache) else{
 bas<-lapply(ths,function(th){cat('constituents',th,'\n');z<-rbindlist(lapply(rr,function(r){if(r-1<2e-4){e<-r-1;return(data.table(y=log(r),logC=log(2)-2+2*log(e),rC=2*r/e,A=(1-log(2))*e,rA=r*(1-log(2))))};cd<-C_log_derivatives(r,th);ad<-action_barrier_derivatives(r,th);data.table(y=log(r),logC=unname(cd['logC']),rC=r*unname(cd['d1']),A=unname(ad['value']),rA=r*unname(ad['d1']))}));d<-1-th
 # Retain the uniform crossover region, then use the asymptotic tail.
 if(d>0){ru<-exp(seq(log(600),log(max(1000,50/d)),length.out=220))[-1];zu<-rbindlist(lapply(ru,function(r){h<-r*2e-4;data.table(y=log(r),logC=large_C(r,th),rC=(large_C(r+h,th)-large_C(r-h,th))/(2*h)*r,A=large_A(r,th),rA=r*large_dA(r,th))}));yy<-seq(log(max(1000,50/d)),max(100,(log(1/min(rhos))+15)/d),length.out=2500)[-1];n<-0:18;zt<-data.table(y=yy,logC=-.5772156649015329-log(d)-yy,rC=-1,A=exp(-d*yy)*vapply(yy,function(y)sum(exp(-n*y)/((n+d)*(n+d+1))),0.),rA=-exp(-d*yy)*vapply(yy,function(y)sum(exp(-n*y)/(n+d+1)),0.));z<-rbind(z,zu,zt)
 }else{ru<-exp(seq(log(600),log(1e8),length.out=1200))[-1];z<-rbind(z,rbindlist(lapply(ru,function(r){h<-r*1e-3;data.table(y=log(r),logC=large_C(r,th),rC=(large_C(r+h,th)-large_C(r-h,th))/(2*h)*r,A=large_A(r,th),rA=r*large_dA(r,th))})))}
 z[,theta:=th];z});saveRDS(bas,cache)
}
all<-list();j<-0
for(th in ths){z<-bas[[which(ths==th)]];iv<-exp(-z$y);lm<-z$y+log1p(-iv)
 for(rv in rhos){Lr<-z$rC+.5/(1-iv)-th/2-z$rA/rv;lb0<-z$logC+.5*(log(rv)+lm-th*z$y-log(2*pi))-z$A/rv
 for(k in ks){lb<-log(k)+lb0;B<-exp(pmin(lb,700));lh<-ifelse(lb< -10,log1p(B/2+B^2/6),ifelse(B<700,log(expm1(pmin(B,700)))-lb,B-lb));val<-Lr+exp(pmin(700,lh-lm));sg<-sign(val);sg[is.na(sg)]<-1;ii<-which(diff(sg)!=0);j<-j+1;all[[j]]<-data.table(rho=rv,theta=th,K=k,nroots=length(ii),first=sg[1],last=tail(sg,1),pattern=paste(c(sg[1],sg[ii+1]),collapse=','))
 }};cat('classified',th,'\n')}
z<-rbindlist(all);fwrite(z,'validation/data/broad_global_structure_scan.csv');print(z[,.N,by=.(nroots,pattern)][order(nroots)]);print(z[first!=1|last!= -1,.N])
