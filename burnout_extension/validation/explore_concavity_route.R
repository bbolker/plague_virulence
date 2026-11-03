source('validation/R/theory.R');source('validation/R/extrema.R');library(data.table)
# Test sufficient condition s_tt<0 on {s>1,s_t<0}; t=log(1-1/R).
rows<-rbindlist(lapply(seq(0,1,.05),function(th)rbindlist(lapply(unique(c(seq(1.01,10,length.out=250),exp(seq(log(10),log(300),length.out=80)))),function(r){
 h<-min(1e-3*(r-1),1e-4*r);cd<-C_log_derivatives(r,th);ad<-action_barrier_derivatives(r,th);dc3<-(C_log_derivatives(r+h,th)['d2']-C_log_derivatives(r-h,th)['d2'])/(2*h)
 c0<-cd['d1']+1/(2*(r-1))-th/(2*r);c1<-cd['d2']-1/(2*(r-1)^2)+th/(2*r^2);c2<-dc3+1/(r-1)^3-th/r^3
 f<-x_final(r);d<-1-r*f;fp<-f*(f-1)/d;de<-1-th;v<-f^de*(de+th*f);vp<-fp*f^(de-1)*(de^2+th*(de+1)*f);dp<--f-r*fp
 a1<-ad['d1'];a2<-ad['d2'];a3<-(th-2)*r^(th-3)/(r-1)-r^(th-2)/(r-1)^2-vp/d+v*dp/d^2
 q<-r*(r-1);qp<-2*r-1
 U<--q*c0;Up<--qp*c0-q*c1;Upp<--2*c0-2*qp*c1-q*c2
 V<-q*a1;Vp<-qp*a1+q*a2;Vpp<-2*a1+2*qp*a2+q*a3
 U1<-q*Up;V1<-q*Vp;U2<-q*(qp*Up+q*Upp);V2<-q*(qp*Vp+q*Vpp)
 lo<-0;hi<-Inf
 # U+V lambda>1; U1+V1 lambda<0
 for(z in list(c(U-1,V),c(-U1,-V1))){if(z[2]>0)lo<-max(lo,-z[1]/z[2])else if(z[2]<0)hi<-min(hi,-z[1]/z[2])else if(z[1]<=0)hi<--1}
 if(hi<=lo)return(NULL)
 bound<-if(V2>0)if(is.infinite(hi))Inf else U2+V2*hi else U2+V2*lo
 if(bound>1e-5){lam<-if(V2>0)if(is.finite(hi))(lo+hi)/2 else max(lo*2,1) else lo*1.01;return(data.table(R0=r,theta=th,lambda=as.numeric(lam),rho=as.numeric(1/lam),s=as.numeric(U+V*lam),s_t=as.numeric(U1+V1*lam),s_tt=as.numeric(U2+V2*lam),upper=as.numeric(bound)))}
 NULL
}))))
fwrite(rows,'validation/data/concavity_sufficient_condition_counterexamples.csv');print(head(rows));cat('counterexample locations',nrow(rows),'\n')
