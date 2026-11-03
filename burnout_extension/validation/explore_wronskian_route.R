source('validation/R/theory.R');source('validation/R/extrema.R');source('validation/R/shape_jets.R');library(data.table)
rr<-unique(c(seq(1.1,12,length.out=100),exp(seq(log(12),log(300),length.out=45))));ths<-seq(0,1,.05);out<-list();j<-0
for(th in ths){for(r in rr){
 z<-shape_derivatives4(r,th);a1<-unname(action_barrier_derivatives(r,th)['d1']);q<-r*(r-1);qp<-2*r-1
 c0<-z['C1']+1/(2*(r-1))-th/(2*r);c1<-z['C2']-1/(2*(r-1)^2)+th/(2*r^2);c2<-z['C3']+1/(r-1)^3-th/r^3;c3<-z['C4']-3/(r-1)^4+3*th/r^4
 U<--q*c0;Up<--qp*c0-q*c1;Upp<--2*c0-2*qp*c1-q*c2;Uppp<--6*c1-3*qp*c2-q*c3
 V<-q*a1;Vp<-qp*a1+q*z['A2'];Vpp<-2*a1+2*qp*z['A2']+q*z['A3'];Vppp<-6*z['A2']+3*qp*z['A3']+q*z['A4']
 U1<-q*Up;V1<-q*Vp;U2<-q*(qp*Up+q*Upp);V2<-q*(qp*Vp+q*Vpp);U3<-q*(qp^2+2*q)*Up+3*q^2*qp*Upp+q^3*Uppp;V3<-q*(qp^2+2*q)*Vp+3*q^2*qp*Vpp+q^3*Vppp
 if(V1<0){j<-j+1;out[[j]]<-data.table(R0=r,theta=th,V1=as.numeric(V1),U1=as.numeric(U1),W=as.numeric(U1*V2-U2*V1),Wt=as.numeric(U1*V3-U3*V1))};next
 lo<-0;hi<-Inf
 for(x in list(c(U-1,V),c(-U1,-V1))){if(x[2]>0)lo<-max(lo,-x[1]/x[2])else if(x[2]<0)hi<-min(hi,-x[1]/x[2])else if(x[1]<=0)hi<--1}
 if(hi<=lo)next
 aa<-V1*V3-V2^2;bb<-U1*V3+V1*U3-2*U2*V2;cc<-U1*U3-U2^2
 vals<-c(lo,if(is.finite(hi))hi else max(lo*2,1));if(aa<0){v<--bb/(2*aa);if(v>lo&&v<hi)vals<-c(vals,v)}
 maxval<-if(is.infinite(hi)&&aa>0)Inf else max(aa*vals^2+bb*vals+cc)
 j<-j+1;out[[j]]<-data.table(R0=r,theta=th,lambda_lo=as.numeric(lo),lambda_hi=as.numeric(hi),quadratic=as.numeric(aa),linear=as.numeric(bb),constant=as.numeric(cc),max_curvature_numerator=as.numeric(maxval))
 };cat('theta',th,'done\n')}
z<-rbindlist(out);fwrite(z,'validation/data/Wronskian_route_checks.csv');print(z[order(Wt)][1:10]);cat('nonpositive Wt',nrow(z[Wt<=0]),'of',nrow(z),'\n')
